
/*********************************************************************************************
 * Copyright (c) <2024>, <Santosh Ansumali@JNCASR>                                                                                                             *
 *  All rights reserved.                                                                                                                                                                          *
 *   Redistribution and use in source and binary forms, with or without modification, are                                                                *
 *   permitted provided that the following conditions are met:                                                                                                           *
 *                                                                                            *
 *    1. Redistributions of source code must retain the above copyright notice, this list of                                                               *
 *       conditions and the following disclaimer.                                                                                                                              *
 *    2. Redistributions in binary form must reproduce the above copyright notice, this list                                                              *
 *       of conditions and the following disclaimer in the documentation and/or other                                                                   *
 *       materials provided with the distribution.                                                                                                                               *
 *    3. Neither the name of the <JNCASR> nor the names of its contributors may be used to                                                         *
 *       endorse or promote products derived from this software without specific prior                                                                  *
 *       written permission.                                                                                                                                                                *
 *                                                                                             *
 *       THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND                           *
 *       ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED                                      *
 *       WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.                     *
 *       IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,                           *
 *       INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,                                  *
 *       BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,                                *
 *       DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF                           *
 *       LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE                                      *
 *       OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED                              *
 *       OF THE POSSIBILITY OF SUCH DAMAGE.                                                                                                                          *
 *                                                                                             *
 *       Suggestions:      ansumali@jncasr.ac.in                                                                                                                  *
 *       Bugs:                 ansumali@jncasr.ac.in                                                                                                                       *
 *       Contributors:                                                                                                                                                                           *
 *       Last Change Date:                                                                                                                                                                  *
 *                                                                                             *
 ************************************************************************************************************************************************/

#ifndef _mesoscale_field_2D_
#define _mesoscale_field_2D_

typedef double real;


#include<iostream>


template<typename dataType, int numField>
class field2D{
    public:
    field2D(int m1, int m2, int numPad){
    
        n1    = m1;
        n2    = m2;
     numDummy = numPad;
        
        
        //Real size of array
        n1R       = m1 + 2 * numPad;
        n2R       = m2 + 2 * numPad;
        numPoints = n1R * n2R;
        
        
        //starting point
        n1Begin    = numPad;
        n2Begin    = numPad;
        
        //end point
        n1End = m1 + numPad - 1;
        n2End = m1 + numPad - 1;
        
        // allocate memory
        //to be written: aligned version of it
    data  = (dataType*) malloc(numField * numPoints * sizeof(dataType));
}

    void makePeriodicX();
    void makePeriodicY();
    
~field2D(){
    // to be written;
}

    int  getIndex(const int i1, const int i2, const int dField)  { return  (i2 * n1R + i1)* numField + dField; }
    dataType  operator()(const int i1, const int i2, const int dField) const     { return data[getIndex(i1,i2,dField)];}
    dataType& operator()(const int i1, const int i2, const int dField)           { return data[getIndex(i1,i2,dField)];}
    dataType  value(const int i1, const int i2, const int dField) const     { return data[getIndex(i1,i2,dField)];}
    dataType& value(const int i1, const int i2, const int dField)           { return data[getIndex(i1,i2,dField)];}
 
 
    int       n1, n2, numDummy, n1R, n2R, numPoints;
    int       n1Begin, n2Begin, n1End, n2End;
    dataType  *data;


};


template<typename dataType, int numField>
void field2D<dataType,numField>::makePeriodicX(){
    
    int wrapD =  n2End- n2Begin+1;
   
    for(int i2 = 0 ; i2 <  n2Begin; i2++){
        for(int i1 = 0 ; i1 < n1R; i1++){
            for(int dv=0;dv< numField;dv++)
                this->value(i1, i2, dv ) = this->value(i1, i2+wrapD, dv );
        }
    }
     
          wrapD =  n2Begin;
       
    for(int i2 =  n2End+1 ; i2 <  n2R; i2++){
        for(int i1 = 0 ; i1 <  n1R; i1++){
            for(int dv=0;dv<numField;dv++)
                this->value(i1, i2, dv ) = this->value(i1,  wrapD, dv );
        }
        wrapD++;
    }

}


template<typename dataType, int numField>
void field2D<dataType,numField>::makePeriodicY(){
    
    int wrapD =  n1End- n1Begin+1;
   
    
        for(int i2 = 0 ; i2 <  n2R; i2++){
            for(int i1 = 0 ; i1 <  n1Begin; i1++){
                for(int dv=0;dv<numField;dv++)
                    this->value(i1, i2, dv ) = this->value(i1+wrapD, i2, dv );
        }
    }
     
          wrapD = n1Begin;
    for(int i2 = 0 ; i2 <  n2R; i2++){
        for(int i1 =  n1End+1 ; i1 < n1R; i1++){
            for(int dv=0;dv<numField;dv++)
                this->value(i1, i2, dv ) = this->value(wrapD, i2, dv );
    }
}
     

}

#endif
