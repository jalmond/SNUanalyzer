// Local include(s):                                                                                                        
#include "SNUError.h"

using namespace std;

/**                                                                                                                         
 * This constructor is used most often in the code. You just specify                                                        
 * a severity to the event, then fill the description using the                                                             
 * std::ostream functionality of the object, finally throw the                                                              
 * exception object.                                                                                                        
 *                                                                                                                          
 * @param severity The action request of the exception                                                                      
 */
SNUError::SNUError( Severity severity ) throw()
  : exception(), ostringstream(), m_severity( severity ) {

}



/**                                                                                                                         
 * This constructor is useful for throwing exceptions with a simple                                                         
 * explanation. For such exceptions you don't have to create an                                                             
 * exception object, but can throw the exception like this:                                                                 
 *                                                                                                                          
 * <code>                                                                                                                   
 *   throw SNUError( "Skip this event", SNUError::SkipEvent );                                                                  
 * </code>                                                                                                                  
 *                                                                                                                          
 * @param description Explanation for the occurance                                                                         
 * @param severity    The action request of the exception                                                                   
 */

SNUError::SNUError( const char* description, Severity severity ) throw()
  : exception(), ostringstream(), m_severity( severity ) {

  this->str( description );

}

/**                                                                                                                         
 * I think this copy constructor is necessary for the correct                                                               
 * handling of the exception objects. It clones the original                                                                
 * object completely.                                                                                                       
 *                                                                                                                          
 * @param parent The object to clone                                                                                        
 */
SNUError::SNUError( const SNUError& parent ) throw()
  : std::basic_ios< SNUError::char_type, SNUError::traits_type >(),
  exception(),
  ostringstream(), m_severity( parent.m_severity ) {

  this->str( parent.str() );

}
/**                                                                                                                         
 * Another "I don't do anything" destructor.                                                                                
 */
SNUError::~SNUError() throw() {

}
void SNUError::SetDescription( const char* description ) throw() {

  this->str( description );
  return;

}

void SNUError::SetSeverity( Severity severity ) throw() {

  m_severity = severity;
  return;

}
const char* SNUError::what() const throw() {

  return this->str().c_str();

}

SNUError::Severity SNUError::request() const throw() {

  return m_severity;

}
