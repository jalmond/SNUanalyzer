// STL include(s):                                                                                                          
#include <map>
#include <cstdlib>

// ROOT include(s):                                                                                                         
#include <TSystem.h>
#include <TError.h>

// Local include(s):  
#include "SNUErrorHandler.h"
#include "SNULogger.h"

/// Local map to translate between ROOT and SNUAnalysis message levels                                                           
static std::map< int, SNUMsgType > msgLevelMap;

/**                                                                                                                         
 * This function is the "SNUAnalysis version" of DefaultErrorHandler defined in the                                              
 * TError.h header. By calling                                                                                              
 *                                                                                                                          
 * <code>                                                                                                                   
 * SetErrorHandler( SNUErrorHandler )                                                                                         
 * </code>                                                                                                                  
 *                                                                                                                          
 * somewhere at the beginning of the application, we can channel all ROOT messages                                          
 * through our own message logging facility.                                                                                
 *                                                                                                                          
 * @param level ROOT message level                                                                                          
 * @param abort Flag telling that the process should abort execution                                                        
 * @param location The source of the message                                                                                
 * @param message The message itself                                                                                        
 */


void SNUErrorHandler( int level, Bool_t abort, const char* location,
                    const char* message ) {

  // Veto some message locations:                                                                                          
  TString tlocation( location );
  if( tlocation.Contains( "NotifyMemory" ) ) {
    return;
  }
  // Create a local logger object:                                                                                         
  SNULogger logger( location );

  // Initialise the helper map the first time the function is called:                                                      
  if( ! msgLevelMap.size() ) {
    msgLevelMap[ kInfo ]     = INFO;
    msgLevelMap[ kWarning ]  = WARNING;
    msgLevelMap[ kError ]    = ERROR;
    msgLevelMap[ kBreak ]    = ERROR;
    msgLevelMap[ kSysError ] = ERROR;
    msgLevelMap[ kFatal ]    = FATAL;
  }


  // Print the message:                                                                                                    
  logger << msgLevelMap[ level ] << message << SNULogger::endmsg;

  // Abort the process if necessary:                                                                                       
  if( abort ) {
    logger << ERROR << "Aborting..." << SNULogger::endmsg;
    if( gSystem ) {
      gSystem->StackTrace();
      gSystem->Abort();
    } else {
      ::abort();
    }
  }

  return;

}
/**                                                                                                                         
 * The following code makes sure that <code>SetErrorHandler(SNUErrorHandler)</code>                                           
 * is called when loading the SNUAnalysisCore library. This way all ROOT messages get                                            
 * printed using SNULogger on the PROOF workers from the moment the SNUAnalysis libraries                                          
 * are loaded. (This is one of the first things that the workers do...)                                                     
 *                                                                                                                          
 * I "stole" the idea for this kind of code from RooFit actually...                                                         
 */
Int_t SetSNUErrorHandler() {

  // Set up SNUAnalysis's error handler:                                                                                        
  SetErrorHandler( SNUErrorHandler );

  // Report this feat:                                                                                                     
  SNULogger logger( "SetSNUErrorHandler" );

  logger << DEBUG << "Redirected ROOT messages to SNUAnalys's logger" << SNULogger::endmsg;

  return 0;

}


// Call the function:                                                                                                       
static Int_t dummy = SetSNUErrorHandler();
