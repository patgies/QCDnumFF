/*
 * ---------------------------------------------------------------------
 *     program testwspace
 * ---------------------------------------------------------------------
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "QCDNUM/QCDNUM.h"
using namespace std;

void K3(double *table, int (&kk)[5]);
int iP3(double *table, int i, int j, int k);
void K2(double *table, int (&kk)[4]);
int iP2(double *table, int i, int j);
void fillhead(double *obj, double *store);
int  sjekhead(double *hd1, double *hd2);

//----------------------------------------------------------------------
int main() {
    int nwords = 573;
    int ntags  = 7;
    int imi[3], ima[3];
    double w[nwords];
    double qq, *q = &qq;
    
    double head[20], hstore[7][20];
    int iws=0, is1=1, it1=2, it2=3, is2=4, it3=5, it4=6;
    double *pws = &hstore[iws][0];
    double *ps1 = &hstore[is1][0];
    double *pt1 = &hstore[it1][0];
    double *pt2 = &hstore[it2][0];
    double *ps2 = &hstore[is2][0];
    double *pt3 = &hstore[it3][0];
    double *pt4 = &hstore[it4][0];
    
    int ndim   =  2;
    imi[0]     =  1;
    ima[0]     =  5;
    imi[1]     =  1;
    ima[1]     =  5;
    
    cout << endl;
    cout << "=============== init workspace ======================" << endl;
    int iaset1 = WSTORE::iws_WsInit(w,nwords,ntags,"Whatta drag");
    cout << endl;
    cout << "=============== add 1st table  ======================" << endl;
    int iatab1 =  WSTORE::iws_WTable(w,imi,ima,ndim);
    cout << endl;
    cout << "=============== add 2nd table  ======================" << endl;
    int iatab2 =  WSTORE::iws_WTable(w,imi,ima,ndim);
    cout << endl;
    cout << "=============== add 2nd tbset  ======================" << endl;
    int iaset2 = WSTORE::iws_NewSet(w);
    cout << endl;
    cout << "=============== add 3rd table  ======================" << endl;
    int iatab3 =  WSTORE::iws_WTable(w,imi,ima,ndim);
    cout << endl;
    cout << "=============== add 4th table  ======================" << endl;
    int iatab4 =  WSTORE::iws_WTable(w,imi,ima,ndim);

    
//  Store headers
    fillhead(w       ,pws);
    fillhead(w+iaset1,ps1);
    fillhead(w+iatab1,pt1);
    fillhead(w+iatab2,pt2);
    fillhead(w+iaset2,ps2);
    fillhead(w+iatab3,pt3);
    fillhead(w+iatab4,pt4);
        
//  Fill table 1
    int index[] = {imi[0],imi[1]};
    int ia = WSTORE::iws_Tpoint(w,iatab1,index,2)-1;
    for(int j=imi[1]; j<=ima[1]; j++) {
      for(int i=imi[0]; i<=ima[0]; i++) {
        ia += 1;
        w[ia] = double(10*i+j);
      }
    }
//  print headers
    WSTORE::sws_WsHead(w,0);
    WSTORE::sws_WsHead(w,iaset1);
    WSTORE::sws_WsHead(w,iatab1);
    WSTORE::sws_WsHead(w,iatab2);
    WSTORE::sws_WsHead(w,iaset2);
    WSTORE::sws_WsHead(w,iatab3);
    WSTORE::sws_WsHead(w,iatab4);
//  print tree
    WSTORE::sws_WsTree(w);
//  Print table 1
    cout << endl << "Contents of table 1" << endl;
    for(int i=imi[0]; i<=ima[0]; i++) {
        for(int j=imi[1]; j<=ima[1]; j++) {
          ia = iP2(w+iatab1,i,j);
          cout << "  " << int(w[ia]) << " ";
        }
      cout << endl;
    }
    
//  Last table and last set address
    int Ltab  = WSTORE::iws_TFskip(w,iatab3);
    int Lset  = WSTORE::iws_SBskip(w,iaset2);
    cout << endl;
    cout << "ia last  table = " << Ltab + iatab3 << endl;
    cout << "ia first tbset = " << Lset + iaset2 << endl;

    cout << endl;
    cout << "Check WsWipe, TsDump and TsRead" << endl;
    cout << "-------------------------------" << endl;
//  Remove last table
    WSTORE::sws_WsWipe(w,iatab4);
    cout << "Wipe  table 4     " << iatab4 << endl;
//  Dump second table-set
    string fname = "../weights/crap.wgt";
    int ierr;
    WSTORE::sws_TsDump(fname,0,w,iaset2,ierr);
    if(ierr!=0) {
        cout << "Error openening or writing file " << fname << endl;
        return 1;
    }
    else {
        cout << "Dump  tableset 2  " << iaset2 << endl;
    }
//  Wipe second table set
    WSTORE::sws_WsWipe(w,iaset2);
    cout << "Wipe  tableset 2  " << iaset2 << endl;
//  Create new (empty) table set
    int iadum = WSTORE::iws_NewSet(w);
    cout << "New   tableset    " << iadum << endl;
//  Read second table set
    iaset2 = WSTORE::iws_TsRead(fname,0,w,ierr);
    if(ierr!=0) {
        cout << "Error openening or reading file " << fname << endl;
        return 1;
    }
    else {
        cout << "Read  tableset 2  " << iaset2 << endl;
    }
//  Restore last table
    iatab4 = WSTORE::iws_WTable(w,imi,ima,ndim);
    cout << "New   table 4     " << iatab4 << endl;
//  Check headers
    cout << endl;
    fillhead(w,head);
    cout << "workspace header OK = " << sjekhead(head,pws) << endl;
    fillhead(w+iaset1,head);
    cout << "set 1 header     OK = " << sjekhead(head,ps1) << endl;
    fillhead(w+iatab1,head);
    cout << "tab 1 header     OK = " << sjekhead(head,pt1) << endl;
    fillhead(w+iatab2,head);
    cout << "tab 2 header     OK = " << sjekhead(head,pt2) << endl;
    fillhead(w+iaset2,head);
    cout << "set 2 header     OK = " << sjekhead(head,ps2) << endl;
    fillhead(w+iatab3,head);
    cout << "tab 3 header     OK = " << sjekhead(head,pt3) << endl;
    fillhead(w+iatab4,head);
    cout << "tab 4 header     OK = " << sjekhead(head,pt4) << endl;
    
    cout << endl;
    cout << "Check TbCopy" << endl;
    cout << "------------" << endl;
    cout << "Contents of table 4 before copy" << endl;
    for(int i=imi[0]; i<=ima[0]; i++) {
        index[0] = i;
        for(int j=imi[1]; j<=ima[1]; j++) {
            index[1] = j;
            ia = WSTORE::iws_Tpoint(w,iatab4,index,2);
            cout << "  " << int(w[ia]) << " ";
        }
        cout << endl;
    }
    cout << "Copy table 1 to table 4" << endl;
    WSTORE::sws_TbCopy(w,iatab1,w,iatab4,1);
    cout << "Contents of table 4 after copy" << endl;
    for(int i=imi[0]; i<=ima[0]; i++) {
        for(int j=imi[1]; j<=ima[1]; j++) {
            ia = iP2(w+iatab4,i,j);
            cout << "  " << int(w[ia]) << " ";
        }
        cout << endl;
    }
    
    cout << endl;
    cout << "Check WClone" << endl;
    cout << "------------" << endl;
//  Wipe second table set
    WSTORE::sws_WsWipe(w,iaset2);
    cout << "Wipe  tableset 2  " << iaset2 << endl;
//  Create new (empty) table set
    iadum = WSTORE::iws_NewSet(w);
    cout << "New   tableset    " << iadum << endl;
//  Clone table set 1
    iaset2 = WSTORE::iws_WClone(w,iaset1,w);
    cout << "Clone set 1 to 2  " << iaset2 << endl;
//  Wipe table 3
    WSTORE::sws_WsWipe(w,iatab3);
    cout << "Wipe  table 3+4   " << iatab3 << endl;
//  Clone table 2 to 3
    iatab3 = WSTORE::iws_WClone(w,iatab2,w);
    cout << "Clone tab 2 to 3  " << iatab3 << endl;
//  Clone table 3 to 4
    iatab4 = WSTORE::iws_WClone(w,iatab3,w);
    cout << "Clone tab 3 to 4  " << iatab4 << endl;
//  Check headers
    cout << endl;
    fillhead(w,head);
    cout << "workspace header OK = " << sjekhead(head,pws) << endl;
    fillhead(w+iaset2,head);
    cout << "set 2 header     OK = " << sjekhead(head,ps2) << endl;
    fillhead(w+iatab3,head);
    cout << "tab 3 header     OK = " << sjekhead(head,pt3) << endl;
    fillhead(w+iatab4,head);
    cout << "tab 4 header     OK = " << sjekhead(head,pt4) << endl;
    
//  Print tree
    WSTORE::sws_WsTree(w);
    
//  Metadata field
    cout << endl;
    cout << "Fingerprint and metadata table1 with address" << iatab1
         << endl;
    cout << "No-can-do because s/r swsGetMeta is Fortran only" << endl;

//  Build-in pointer function
    cout << endl;
    cout << "Test iws_Tpoint function on table1" << endl;
    index[0] =  imi[0];
    index[1] =  imi[1];
    ia     = WSTORE::iws_Tpoint(w,iatab1,index,2);
    int ib = WSTORE::iws_BeginTbody(w,iatab1);
    cout << "begin table = " << ia << " " << ib << endl;
    index[0] =  ima[0];
    index[1] =  ima[1];
    ia = WSTORE::iws_Tpoint(w,iatab1,index,2);
    ib = WSTORE::iws_EndTbody(w,iatab1);
    cout << "end   table = " << ia << " " << ib << endl;
    
    int icWsp, icSet, icTab;
    WSTORE::sws_WsMark(icWsp,icSet,icTab);
    
    cout << endl;
    cout << "tsize  = " << WSTORE::iws_TbSize(imi,ima,ndim) << endl;
    cout << "hsize  = " << WSTORE::iws_HdSize() << endl;
    cout << "Wspace = " << icWsp << endl;
    cout << "TbSet  = " << icSet << endl;
    cout << "Table  = " << icTab << endl;
    cout << "Yesno  = " << WSTORE::iws_IsaWorkspace(head) << endl;
    cout << "Wsize  = " << WSTORE::iws_SizeOfW(w) << endl;
    cout << "Wused  = " << WSTORE::iws_WordsUsed(w) << endl;
    cout << "Nhead  = " << WSTORE::iws_Nheader(w) << endl;
    cout << "Ntags  = " << WSTORE::iws_Ntags(w) << endl;
    cout << "Hskip  = " << WSTORE::iws_HeadSkip(w) << endl;
    cout << "Drain  = " << WSTORE::iws_IaDrain(w) << endl;
    cout << "INull  = " << WSTORE::iws_IaNull(w) << endl;
    int ityp0 = WSTORE::iws_ObjectType(w,66);
    int ityp1 = WSTORE::iws_ObjectType(w,0);
    int ityp2 = WSTORE::iws_ObjectType(w,iaset1);
    int ityp3 = WSTORE::iws_ObjectType(w,iatab2);
    cout << "Itype  = " << ityp0 << " " << ityp1 << " "
                        << ityp2 << " " << ityp3 << endl;
    cout << "WSsize = " << WSTORE::iws_ObjectSize(w,0) << endl;
    cout << "TSsize = " << WSTORE::iws_ObjectSize(w,iaset1) << endl;
    cout << "TBsize = " << WSTORE::iws_ObjectSize(w,iatab1) << endl;
    cout << "Nsets  = " << WSTORE::iws_Nobjects(w,0) << endl;
    cout << "Ntabs  = " << WSTORE::iws_ObjectNumber(w,iaset1) << " "
                        << WSTORE::iws_Nobjects(w,iaset1) << endl;
    cout << "Ntabs  = " << WSTORE::iws_ObjectNumber(w,iaset2) << " "
                        << WSTORE::iws_Nobjects(w,iaset2) << endl;
    cout << "Wfprnt = " << WSTORE::iws_FingerPrint(w,0) << endl;
    cout << "S1fp   = " << WSTORE::iws_FingerPrint(w,iaset1) << endl;
    cout << "S2fp   = " << WSTORE::iws_FingerPrint(w,iaset2) << endl;
    cout << "T1fp   = " << WSTORE::iws_FingerPrint(w,iatab1) << endl;
    cout << "T2fp   = " << WSTORE::iws_FingerPrint(w,iatab2) << endl;
    cout << "Wtag   = " << WSTORE::iws_IaFirstTag(w,0) << endl;
    cout << "Ndim   = " << WSTORE::iws_TableDim(w,iatab2) << endl;
    cout << "KARR1  = " << WSTORE::iws_IaKARRAY(w,iatab1) << endl;
    cout << "IMIN1  = " << WSTORE::iws_IaIMIN(w,iatab1) << endl;
    cout << "IMAX1  = " << WSTORE::iws_IaIMAX(w,iatab1) << endl;
    cout << "Body1  = " << WSTORE::iws_BeginTbody(w,iatab1) << endl;
    cout << "Bend1  = " << WSTORE::iws_EndTbody(w,iatab1) << endl;
    
    return 0;
}

//----------------------------------------------------------------------

void fillhead(double *obj, double *store) {
    int n = WSTORE::iws_HdSize();
    MBUTIL::smb_vcopy(obj,store,n);
}

int  sjekhead(double *hd1, double *hd2) {
    int n = WSTORE::iws_HdSize();
    return MBUTIL::lmb_vcomp(hd1,hd2,n,-1e-9);
}

void K3(double *table, int (&kk)[5]) {
    int ia    = int(*(table+1));
    double *w = table-ia;
    int iak   = WSTORE::iws_IaKARRAY(w,ia);
    *(kk)     = WSTORE::iws_FingerPrint(w,ia);
    *(kk+1)   = int(*(w+iak));
    *(kk+2)   = int(*(w+iak+1));
    *(kk+3)   = int(*(w+iak+2));
    *(kk+4)   = int(*(w+iak+3));
}

int iP3(double *table, int i, int j, int k) {
    static int kk[5];
    int ia    = int(*(table+1));
    double *w = table-ia;
    int ifp   = WSTORE::iws_FingerPrint(w,ia);
    if( kk[0] != ifp ) { K3( table, kk ); }
    int ip    = kk[1]+kk[2]*i+kk[3]*j+kk[4]*k;
    return ia + ip;
}
                              
void K2(double *table, int (&kk)[4]) {
    int ia    = int(*(table+1));
    double *w = table-ia;
    int iak   = WSTORE::iws_IaKARRAY(w,ia);
    *(kk)     = WSTORE::iws_FingerPrint(w,ia);
    *(kk+1)   = int(*(w+iak));
    *(kk+2)   = int(*(w+iak+1));
    *(kk+3)   = int(*(w+iak+2));
}

int iP2(double *table, int i, int j) {
    static int kk[4];
    int ia    = int(*(table+1));
    double *w = table-ia;
    int ifp   = WSTORE::iws_FingerPrint(w,ia);
    if( kk[0] != ifp ) { K2( table, kk ); }
    int ip    = kk[1]+kk[2]*i+kk[3]*j;
    return ia + ip;
}
                              
