#include <stdio.h>
#include "CoolPropC.h"
//   std::string str_out(Output,getLength(Output));
//   std::string str_ref(Ref,getLength(Ref));
// //   printf("String for Output: %s\n",str_out.c_str());
// //   printf("String for Ref   : %s\n",str_ref.c_str());
//   
// //   //printf("double for res : %4.1f\n",Props(str_out.c_str(),(char *)"T",&T,(char *)"Q",&Q,str_ref.c_str()));
// //   //printf("double for res : %4.1f\n",Props(Output, *Name1, *Prop1, *Name2, *Prop2, Ref));
// //   double res = Props(str_out, *Name1, *Prop1, *Name2, *Prop2, str_ref);
// //   printf("double for res : %4.1f\n",res);
// //   //*outVal = res;
//   return Props(str_out.c_str(), *Name1, *Prop1, *Name2, *Prop2, (char *) str_ref.c_str());
// }
int main(int argc, char *argv[]) 
{
  double T =  285;
  double Q =    0;
  double D = 1250;
  printf("double for props1    : %4.1f\n",              props1_((char *)"R134a" , (char *)"Tcrit"));
  printf("double for props     : %4.1f\n",               props_((char *)"S"     , (char *)"T"    , &T, (char *)"Q",    &Q, (char *)"R134a"));
  printf("double for deriv     : %4.1f\n",          derivterms_((char *)"dpdrho", &T             , &D, (char *)"R134a"));
  printf("int for reference    : %d   \n",set_reference_states_((char *)"R134a" , (char *)"IIR"));
  printf("double for props     : %4.1f\n",               props_((char *)"S"     , (char *)"T"    , &T, (char *)"Q",    &Q, (char *)"R134a"));
  printf("int for reference    : %d   \n",set_reference_states_((char *)"R134a" , (char *)"ASHRAE"));
  printf("double for props     : %4.1f\n",               props_((char *)"S"     , (char *)"T"    , &T, (char *)"Q",    &Q, (char *)"R134a"));
  printf("bool for enable ttse : %d   \n",     enable_ttse_lut_((char *)"R134a"));
  printf("bool for check ttse  : %d   \n",  isenabled_ttse_lut_((char *)"R134a"));
  printf("double for props     : %4.1f\n",               props_((char *)"S"     , (char *)"T"    , &T, (char *)"Q",    &Q, (char *)"R134a"));
  printf("int for set ttse mode: %d   \n",       set_ttse_mode_((char *)"R134a" , (char *)"bicubic"));
  printf("double for props     : %4.1f\n",               props_((char *)"S"     , (char *)"T"    , &T, (char *)"Q",    &Q, (char *)"R134a"));  
  printf("bool for disable ttse: %d   \n",    disable_ttse_lut_((char *)"R134a"));
  printf("double for props     : %4.1f\n",               props_((char *)"S"     , (char *)"T"    , &T, (char *)"Q",    &Q, (char *)"R134a"));
  return(1);
}