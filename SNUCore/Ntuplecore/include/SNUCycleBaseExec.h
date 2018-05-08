#ifndef SNUCycleBaseExec_H
#define SNUCycleBaseExec_H

// Local include(s):
#include "SNUCycleBaseBase.h"
#include "SKTreeFiller.h"

class SNUCycleBaseExec :   public virtual SNUCycleBaseBase , public SKTreeFiller{

 public:
  /// Default constructor  
  SNUCycleBaseExec();

  ///////////////////////////////////////////////////////////////////////////
  //                                                                       //
  //   The following are the functions to be implemented in the derived    //
  //   classes.                                                            //
  //                                                                       //
  ///////////////////////////////////////////////////////////////////////////

  /// Initialisation called at the beginning of a full cycle
  /**
   * Analysis-wide configurations, like the setup of some reconstruction
   * algorithm based on properties 
   */

  virtual void BeginCycle() throw( SNUError );

  /**
   * Called before the  event. Gets the weight from the configured job
   **/
  virtual void BeginEvent()throw( SNUError );
  
  /// Function called for every event
  /**
   * This is the function where the main analysis should be done. By the
   * time it is called, all the input variables are filled with the
   * contents of the actual event.
   */
  virtual void ExecuteEvents() throw( SNUError );

  /**  
   *  Called at end of cycle
   *
   **/
  virtual void EndCycle()throw( SNUError );
  virtual void EndEvent()throw( SNUError );
  
  /**
   *   Interates through the ntuple
   *
   **/
  virtual void SetUpEvent(Long64_t jevent, float ev_weight)throw( SNUError );

  /**
   *  call before each event. Specific function to clear output vectors
   *
   **/

  virtual void ClearOutputVectors() throw( SNUError );

  /** 
   *  call al end of cycle: Will write all histograms to output file
   **/

  virtual void WriteHistograms() throw( SNUError );
  
  ClassDef(SNUCycleBaseExec, 0 );
};

#endif
