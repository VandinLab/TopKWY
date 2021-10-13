// patterngraph.cpp
// Siegfried Nijssen, snijssen@liacs.nl, jan 2004.
#include "patterngraph.h"
#include "graphstate.h"

void PatternGraph::init ( vector<CloseLegPtr> &closelegssource, int legindex ) {
  #ifdef deepdebug
  cout << "call to PatternGraph::init" << endl;
  #endif

  closetuples.push_back ( closelegssource[legindex]->tuple );
  frequency = closelegssource[legindex]->occurrences.frequency;
  this->current_info.frequency = frequency;
  this->current_info.max_as = -1;
  this->current_info.min_as = -1;
  this->current_info.as = -1;

  // ADDED
  graphstate.closetuples = &closetuples;

  this->closelegssource = &closelegssource;
  this->legindex = legindex;
  #ifdef deepdebug
  cout << "call to PatternGraph::init done" << endl;
  #endif

}


extern bool f;

void PatternGraph::expand () {
  #ifdef deepdebug
  cout << "call to PatternGraph::expand" << endl;
  #endif

  // ADDED
  int id = graphstate.isnormal ();
  if ( id == 0 ) {
    for ( int k = legindex + 1; k < closelegssource->size (); k++ ) {
      if ( (*closelegssource)[k]->copy ) {
        CloseLegOccurrencesPtr closelegoccurrencesptr =
          join ( (*closelegssource)[legindex]->occurrences, (*closelegssource)[k]->occurrences );
        if ( closelegoccurrencesptr ) {
          CloseLegPtr closelegptr = new CloseLeg;
          closelegptr->tuple = (*closelegssource)[k]->tuple;
          swap ( *closelegoccurrencesptr, closelegptr->occurrences );
          closelegs.push_back ( closelegptr );
        }
      }
    }
    // OUTPUT(frequency)
    #ifdef deepdebug
    cout << "   processing CL graph of support " << (*closelegssource)[legindex]->occurrences.frequency << endl;
    cout << "   support of parent " << this->parent_info.frequency << endl; // correct
    #endif
    bounds_info current_info_temp;
    current_info_temp.frequency = (*closelegssource)[legindex]->occurrences.frequency;
    checkTestableCl((*closelegssource)[legindex]->occurrences.elements, (*closelegssource)[legindex]->occurrences.frequency, this->parent_info, current_info_temp);
    #ifdef deepdebug
    cout << "   done processing CL graph of support " << (*closelegssource)[legindex]->occurrences.frequency << endl;
    #endif
    int addsize = statistics.patternsize + graphstate.edgessize - graphstate.nodes.size ();
    if ( addsize >= statistics.frequenttreenumbers.size () ) {
      statistics.frequenttreenumbers.resize ( addsize + 1, 0 );
      statistics.frequentpathnumbers.resize ( addsize + 1, 0 );
      statistics.frequentgraphnumbers.resize ( addsize + 1, 0 );
    }
    statistics.frequentgraphnumbers[addsize]++;

    for ( int k = closelegs.size () - 1; k >= 0; k-- ) {
      graphstate.insertEdge ( closelegs[k]->tuple.from, closelegs[k]->tuple.to, closelegs[k]->tuple.label );
      PatternGraph patterngraph ( *this, k );
      patterngraph.parent_info = this->current_info;
      //patterngraph.parent_frequency = this->frequency;
      patterngraph.expand ();
      graphstate.deleteEdge ( closelegs[k]->tuple.from, closelegs[k]->tuple.to );
    }
  }
  else
    if ( id == 2 ) {
      (*closelegssource)[legindex]->copy = false; // should not be added to any later tree
    }
    #ifdef deepdebug
    cout << "call to PatternGraph::expand done" << endl;
    #endif
}

PatternGraph::~PatternGraph () {
  for ( int i = 0; i < closelegs.size (); i++ )
    delete closelegs[i];
}
