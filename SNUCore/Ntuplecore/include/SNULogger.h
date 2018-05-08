#ifndef SNULogger_H
#define SNULogger_H

// STL include(s):
#include <string>
#include <sstream>

// Local include(s):
#include "SNUMsgType.h"
#include "SNULogWriter.h"

// Forward declaration(s):
class TObject;
class SNULogWriter;

class SNULogger : public std::ostringstream {
  
 public:
  /// Constructor with pointer to the parent object
  SNULogger( const TObject* source );
  /// Constructor with a name of the parent object
  SNULogger( const std::string& source );
  /// Copy constructor
  SNULogger( const SNULogger& parent );
  /// Default destructor
  virtual ~SNULogger();

  void SetSource( const TObject* source );
  void SetSource( const std::string& source );

  /// Copy operator
  SNULogger& operator= ( const SNULogger& parent );

  /// Stream modifier to send a message
  static SNULogger& endmsg( SNULogger& logger );

  /// Operator accepting SNULogger stream modifiers
  SNULogger& operator<< ( SNULogger& ( *_f )( SNULogger& ) );
  /// Operator accepting std::ostream stream modifiers
  SNULogger& operator<< ( std::ostream& ( *_f )( std::ostream& ) );
  /// Operator accepting std::ios stream modifiers
  SNULogger& operator<< ( std::ios& ( *_f )( std::ios& ) );

  /// Operator accepting message type setting
  SNULogger& operator<< ( SNUMsgType type );

  /// Operator accepting basically any kind of argument
  /**
   * SNULogger was designed to give all the features that std::ostream
   * objects usually provide. This operator handles all kinds of
   * arguments and passes it on to the std::ostringstream base class.
   */
  template < class T > SNULogger& operator<< ( T arg ) {
    if( m_activeType >= m_logWriter->GetMinType() ) {
      ( * ( std::ostringstream* ) this ) << arg;
    }
    return *this;
  }

  /// Old style message sender function
  void Send( SNUMsgType type, const std::string& message ) const;

 private:
  void Send();

  const TObject* m_objSource;
  std::string    m_strSource;


  SNULogWriter*    m_logWriter;
  SNUMsgType       m_activeType;

}; // class SNULogger

//////////////////////////////////////////////////////////////////////
//                                                                  //
//   To speed up the code a bit, the following operators are        //
//   declared 'inline'.                                             //
//                                                                  //
//////////////////////////////////////////////////////////////////////

/**
 * This operator handles all stream modifiers that have been written
 * to work on SNULogger objects specifically. Right now there is basically
 * only the SNULogger::endmsg stream modifier that is such.
 */
inline SNULogger& SNULogger::operator<< ( SNULogger& ( *_f )( SNULogger& ) ) {

  return ( _f )( *this );
}

/**
 * This operator handles all stream modifiers that have been written
 * to work on std::ostream objects. Most of the message formatting
 * modifiers are such.
 */
inline SNULogger& SNULogger::operator<< ( std::ostream& ( *_f )( std::ostream& ) ) {

  if( m_activeType >= m_logWriter->GetMinType() ) {
    ( _f )( *this );
  }
  return *this;
}


/**
 * This operator handles all stream modifiers that have been written
 * to work on std::ios objects. I have to admit I don't remember exactly
 * which operators these are, but some formatting operations need this.
 */

inline SNULogger& SNULogger::operator<< ( std::ios& ( *_f )( std::ios& ) ) {

  if( m_activeType >= m_logWriter->GetMinType() ) {
    ( _f )( *this );
  }
  return *this;
}

/**
 * Messages have a type, defined by the SNUMsgType enumeration. This operator
 * allows the user to write intuitive message lines in the code like this:
 *
 * <code>
 *   logger << INFO << "This is an info message" << SNULogger::endmsg;
 * </code>
 */
inline SNULogger& SNULogger::operator<< ( SNUMsgType type ) {

  m_activeType = type;
  return *this;
}

// This is a GCC extension for getting the name of the current function.
#if defined( __GNUC__ )
#   define SNULOGGER_FNAME __PRETTY_FUNCTION__
#else
#   define SNULOGGER_FNAME ""
#endif

/// Common prefix for the non-usual messages
/**
 * The idea is that a regular user usually only wants to see DEBUG, INFO
 * and some WARNING messages. So those should be reasonably short. On the other
 * hand serious warnings (ERROR, FATAL) or VERBOSE messages should be as precise
 * as possible.
 *
 * So I stole the idea from Athena (what a surprise...) to have a few macros which
 * produce messages with a common formatting. This macro provides the prefix for
 * all the messages.
 */
#define SNULOGGER_REPORT_PREFIX \
   __FILE__ << ":" << __LINE__ << " (" << SNULOGGER_FNAME << "): "

/// Convenience macro for reporting VERBOSE messages in the code
/**
 * This macro is very similar to the REPORT_MESSAGE macros of Athena. It prints
 * a nicely formatted output that specifies both the exact function name where
 * the message was printed, and also the filename:line combination. It can be used
 * like a regular function inside cycles:
 *
 * <code>
 *   REPORT_VERBOSE( "This is a verbose message with a number: " << number );
 * </code>
 */
#define REPORT_VERBOSE( MESSAGE ) \
  m_logger << VERBOSE << SNULOGGER_REPORT_PREFIX << MESSAGE << SNULogger::endmsg

/// Convenience macro for reporting ERROR messages in the code
/**
 * This macro is very similar to the REPORT_MESSAGE macros of Athena. It prints
 * a nicely formatted output that specifies both the exact function name where
 * the message was printed, and also the filename:line combination. It can be used
 * like a regular function inside cycles:
 *
 * <code>
 *   REPORT_ERROR( "A serious error message" );
 * </code>
 */
#define REPORT_ERROR( MESSAGE ) \
   m_logger << ERROR << SNULOGGER_REPORT_PREFIX << MESSAGE << SNULogger::endmsg

/// Convenience macro for reporting FATAL messages in the code

/**
* This macro is very similar to the REPORT_MESSAGE macros of Athena. It prints
* a nicely formatted output that specifies both the exact function name where
* the message was printed, and also the filename:line combination. It can be used
* like a regular function inside cycles:
*
* <code>
*   REPORT_FATAL( "A very serious error message" );
* </code>
*/
#define REPORT_FATAL( MESSAGE ) \
  m_logger << FATAL << SNULOGGER_REPORT_PREFIX << MESSAGE << SNULogger::endmsg

#endif // SNULogger_H
