#ifndef SNULogWriter_H
#define SNULogWriter_H


// STL include(s):
#include <string>
#include <map>

// Local include(s):
#include "SNUMsgType.h"


/**
 *   @short Message writing class
 *
 *          Singleton class for actually writing the formatted
 *          messages to the console.
 *
 *          Right now it only writes messages to the terminal, but
 *          one possibility would be to write messages to a file
 *          for batch running later on. (Just an idea...)
 *
 *     @see SNULogger
 * @version $Revision: 1 $
 */

class SNULogWriter {

 public:
  /// Function for accessing the single object
  static SNULogWriter* Instance();
  /// Default destructor
  ~SNULogWriter();

  /// Function writing a message to the output
  void Write( SNUMsgType type, const std::string& line ) const;

  /// Set the message type above which messages are printed
  void SetMinType( SNUMsgType type );
  /// Get the message type above which messages are printed
  SNUMsgType GetMinType() const;

 protected:
  /// Protected default constructor
  SNULogWriter();

 private:
  static SNULogWriter* m_instance;

  std::map< SNUMsgType, std::string > m_typeMap;
  std::map< SNUMsgType, std::string > m_colorMap;
  SNUMsgType                          m_minType;

}; // class SNULogWriter

#endif //SNULogWriter_H
