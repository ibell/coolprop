#include <string.h>
#include <stdio.h>
#include <CoolProp.h>

#ifdef __cplusplus
extern "C" {
#endif
double props_(char *Output, char *Name1, double *Prop1, char *Name2, double *Prop2, char *Ref);
#ifdef __cplusplus
}
#endif

int getLength(char *string)
{
   //chars[theLength--] = '\0';  // NULL terminate the string
   int len = strlen(string);
   while (string[len-1] == ' ' || string[len-1] == '\0') --len;
   return len;
}

//int props_(char *Output, char *Name1, double Prop1, char *Name2, double Prop2, char *Ref, double *outVal, int outLen, int name1Len, int name2Len, int refLen)
//int props_(char Output[fString], char Name1[1])
double props_(char *Output, char *Name1, double *Prop1, char *Name2, double *Prop2, char *Ref)//, int outLen, int name1Len, int name2Len, int refLen)
{
  
//   printf("char for Output  : %s\n",Output);
//   printf("char for Name1   : %s\n",Name1);
//   printf("char for Name2   : %s\n",Name2);
//   printf("char for Ref     : %s\n",Ref);
//   printf("Length for Output: %d\n",strlen(Output));
//   printf("Length for Name1 : %d\n",strlen(Name1));
//   printf("Length for Name2 : %d\n",strlen(Name2));
//   printf("Length for Ref   : %d\n",strlen(Ref));
//   printf("Length for Output: %d\n",getLength(Output));
//   printf("Length for Name1 : %d\n",getLength(Name1));
//   printf("Length for Name2 : %d\n",getLength(Name2));
//   printf("Length for Ref   : %d\n",getLength(Ref));
//   printf("double for Prop1 : %4.1f\n",*Prop1);
//   printf("double for Prop2 : %4.1f\n",*Prop2);
// // //   printf("Length for Output: %d\n",outLen);
// // //   printf("Length for Name1 : %d\n",name1Len);
// // //   printf("Length for Name2 : %d\n",name2Len);
// // //   printf("Length for Ref   : %d\n",refLen);
// //   printf("String for Output: %.*s\n",outLen,Output);
// //   printf("String for Name1 : %.*s\n",name1Len,Name1);
// //   printf("String for Name2 : %.*s\n",name2Len,Name2);
// //   printf("String for Ref   : %.*s\n",refLen,Ref);
  
  // Process the inputs
  std::string str_out(Output,getLength(Output));
  std::string str_ref(Ref,getLength(Ref));
//   printf("String for Output: %s\n",str_out.c_str());
//   printf("String for Ref   : %s\n",str_ref.c_str());
  
//   //printf("double for res : %4.1f\n",Props(str_out.c_str(),'T',285.0,'Q',0.0,str_ref.c_str()));
//   //printf("double for res : %4.1f\n",Props(Output, *Name1, *Prop1, *Name2, *Prop2, Ref));
//   double res = Props(str_out, *Name1, *Prop1, *Name2, *Prop2, str_ref);
//   printf("double for res : %4.1f\n",res);
//   //*outVal = res;
  return Props(str_out.c_str(), *Name1, *Prop1, *Name2, *Prop2, (char *) str_ref.c_str());
}
