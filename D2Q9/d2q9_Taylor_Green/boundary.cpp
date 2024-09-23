/*********************************************************************************************
 * Copyright (c) 2024, Santosh Ansumali@JNCASR
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without modification, 
 * are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice, 
 *    this list of conditions and the following disclaimer in the documentation 
 *    and/or other materials provided with the distribution.
 * 3. Neither the name of JNCASR nor the names of its contributors may be used to 
 *    endorse or promote products derived from this software without specific 
 *    prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
 * POSSIBILITY OF SUCH DAMAGE.
 * 
 * Suggestions: ansumali@jncasr.ac.in
 * Bugs: ansumali@jncasr.ac.in
 * Contributors:
 * Last Change Date:
 *********************************************************************************************/

template <typename dataType>
void periodicX(field2D<dataType, 9> &lbmGrid, lbmD2Q9<dataType> &d2q9Model) {
  for (int i2 = 0; i2 < lbmGrid.n2R; i2++) {
    for (int dv = 0; dv < 9; dv++) {
      lbmGrid(0, i2, dv) = lbmGrid(lbmGrid.n1End, i2, dv);
      lbmGrid(lbmGrid.n1End + 1, i2, dv) = lbmGrid(1, i2, dv);
    }
  }
}

template <typename dataType>
void periodicY(field2D<dataType, 9> &lbmGrid, lbmD2Q9<dataType> &d2q9Model) {
  for (int i1 = 0; i1 < lbmGrid.n1R; i1++) {
    for (int dv = 0; dv < 9; dv++) {
      lbmGrid(i1, 0, dv) = lbmGrid(i1, lbmGrid.n2End, dv);
      lbmGrid(i1, lbmGrid.n2End + 1, dv) = lbmGrid(i1, 1, dv);
    }
  }
}