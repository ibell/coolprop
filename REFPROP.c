#if defined (__ISWINDOWS__) //Check if it is a windows machine, if not, hide this function
double REFPROP(char Output,char Name1, double Prop1, char Name2, double Prop2, char * Ref)
{
	int j;
	long i,ierr=0;
	char hf[refpropcharlength*ncmax], hrf[lengthofreference+1],
	herr[errormessagelength+1],hfmix[refpropcharlength+1];
	
	double x[ncmax],xliq[ncmax],xvap[ncmax];
	char RefString[255];
	double T,p=0,d,dl,dv,q,e,h,s,cv,cp,w,MW,hl,hv,sl,sv,ul,uv,pl,pv,eta,tcx,Q,Tcrit,pcrit,dcrit,rho;

	// First create a pointer to an instance of the library
	// Then have windows load the library.
		HINSTANCE RefpropdllInstance;
		RefpropdllInstance = LoadLibrary("c:/Program Files/REFPROP/refprop.dll");

	// Then get pointers into the dll to the actual functions.
		ABFL1dll = (fp_ABFL1dllTYPE) GetProcAddress(RefpropdllInstance,"ABFL1dll");
		ABFL2dll = (fp_ABFL2dllTYPE) GetProcAddress(RefpropdllInstance,"ABFL2dll");
		ACTVYdll = (fp_ACTVYdllTYPE) GetProcAddress(RefpropdllInstance,"ACTVYdll");
		AGdll = (fp_AGdllTYPE) GetProcAddress(RefpropdllInstance,"AGdll");
		CCRITdll = (fp_CCRITdllTYPE) GetProcAddress(RefpropdllInstance,"CCRITdll");
		CP0dll = (fp_CP0dllTYPE) GetProcAddress(RefpropdllInstance,"CP0dll");
		CRITPdll = (fp_CRITPdllTYPE) GetProcAddress(RefpropdllInstance,"CRITPdll");
		CSATKdll = (fp_CSATKdllTYPE) GetProcAddress(RefpropdllInstance,"CSATKdll");
		CV2PKdll = (fp_CV2PKdllTYPE) GetProcAddress(RefpropdllInstance,"CV2PKdll");
		CVCPKdll = (fp_CVCPKdllTYPE) GetProcAddress(RefpropdllInstance,"CVCPKdll");
		CVCPdll = (fp_CVCPdllTYPE) GetProcAddress(RefpropdllInstance,"CVCPdll"); 
		DBDTdll = (fp_DBDTdllTYPE) GetProcAddress(RefpropdllInstance,"DBDTdll");
		DBFL1dll = (fp_DBFL1dllTYPE) GetProcAddress(RefpropdllInstance,"DBFL1dll");
		DBFL2dll = (fp_DBFL2dllTYPE) GetProcAddress(RefpropdllInstance,"DBFL2dll");
		DDDPdll = (fp_DDDPdllTYPE) GetProcAddress(RefpropdllInstance,"DDDPdll");
		DDDTdll = (fp_DDDTdllTYPE) GetProcAddress(RefpropdllInstance,"DDDTdll");
		DEFLSHdll = (fp_DEFLSHdllTYPE) GetProcAddress(RefpropdllInstance,"DEFLSHdll");
		DHD1dll = (fp_DHD1dllTYPE) GetProcAddress(RefpropdllInstance,"DHD1dll");
		DHFLSHdll = (fp_DHFLSHdllTYPE) GetProcAddress(RefpropdllInstance,"DHFLSHdll");
		DIELECdll = (fp_DIELECdllTYPE) GetProcAddress(RefpropdllInstance,"DIELECdll");
		DOTFILLdll = (fp_DOTFILLdllTYPE) GetProcAddress(RefpropdllInstance,"DOTFILLdll");
		DPDD2dll = (fp_DPDD2dllTYPE) GetProcAddress(RefpropdllInstance,"DPDD2dll");
		DPDDKdll = (fp_DPDDKdllTYPE) GetProcAddress(RefpropdllInstance,"DPDDKdll");
		DPDDdll = (fp_DPDDdllTYPE) GetProcAddress(RefpropdllInstance,"DPDDdll");
		DPDTKdll = (fp_DPDTKdllTYPE) GetProcAddress(RefpropdllInstance,"DPDTKdll");
		DPDTdll = (fp_DPDTdllTYPE) GetProcAddress(RefpropdllInstance,"DPDTdll");
		DPTSATKdll = (fp_DPTSATKdllTYPE) GetProcAddress(RefpropdllInstance,"DPTSATKdll");
		DSFLSHdll = (fp_DSFLSHdllTYPE) GetProcAddress(RefpropdllInstance,"DSFLSHdll");
		ENTHALdll = (fp_ENTHALdllTYPE) GetProcAddress(RefpropdllInstance,"ENTHALdll"); //**
		ENTROdll = (fp_ENTROdllTYPE) GetProcAddress(RefpropdllInstance,"ENTROdll");
		ESFLSHdll = (fp_ESFLSHdllTYPE) GetProcAddress(RefpropdllInstance,"ESFLSHdll");
		FGCTYdll = (fp_FGCTYdllTYPE) GetProcAddress(RefpropdllInstance,"FGCTYdll");
		FPVdll = (fp_FPVdllTYPE) GetProcAddress(RefpropdllInstance,"FPVdll");
		GERG04dll = (fp_GERG04dllTYPE) GetProcAddress(RefpropdllInstance,"GERG04dll");
		GETFIJdll = (fp_GETFIJdllTYPE) GetProcAddress(RefpropdllInstance,"GETFIJdll");
		GETKTVdll = (fp_GETKTVdllTYPE) GetProcAddress(RefpropdllInstance,"GETKTVdll");
		GIBBSdll = (fp_GIBBSdllTYPE) GetProcAddress(RefpropdllInstance,"GIBBSdll");
		HSFLSHdll = (fp_HSFLSHdllTYPE) GetProcAddress(RefpropdllInstance,"HSFLSHdll");
		INFOdll = (fp_INFOdllTYPE) GetProcAddress(RefpropdllInstance,"INFOdll");
		LIMITKdll = (fp_LIMITKdllTYPE) GetProcAddress(RefpropdllInstance,"LIMITKdll");
		LIMITSdll = (fp_LIMITSdllTYPE) GetProcAddress(RefpropdllInstance,"LIMITSdll");
		LIMITXdll = (fp_LIMITXdllTYPE) GetProcAddress(RefpropdllInstance,"LIMITXdll");
		MELTPdll = (fp_MELTPdllTYPE) GetProcAddress(RefpropdllInstance,"MELTPdll");
		MELTTdll = (fp_MELTTdllTYPE) GetProcAddress(RefpropdllInstance,"MELTTdll");
		MLTH2Odll = (fp_MLTH2OdllTYPE) GetProcAddress(RefpropdllInstance,"MLTH2Odll");
		NAMEdll = (fp_NAMEdllTYPE) GetProcAddress(RefpropdllInstance,"NAMEdll");
		PDFL1dll = (fp_PDFL1dllTYPE) GetProcAddress(RefpropdllInstance,"PDFL1dll");
		PDFLSHdll = (fp_PDFLSHdllTYPE) GetProcAddress(RefpropdllInstance,"PDFLSHdll");
		PEFLSHdll = (fp_PEFLSHdllTYPE) GetProcAddress(RefpropdllInstance,"PEFLSHdll");
		PHFL1dll = (fp_PHFL1dllTYPE) GetProcAddress(RefpropdllInstance,"PHFL1dll");
		PHFLSHdll = (fp_PHFLSHdllTYPE) GetProcAddress(RefpropdllInstance,"PHFLSHdll");
		PQFLSHdll = (fp_PQFLSHdllTYPE) GetProcAddress(RefpropdllInstance,"PQFLSHdll");
		PREOSdll = (fp_PREOSdllTYPE) GetProcAddress(RefpropdllInstance,"PREOSdll");
		PRESSdll = (fp_PRESSdllTYPE) GetProcAddress(RefpropdllInstance,"PRESSdll");
		PSFL1dll = (fp_PSFL1dllTYPE) GetProcAddress(RefpropdllInstance,"PSFL1dll");
		PSFLSHdll = (fp_PSFLSHdllTYPE) GetProcAddress(RefpropdllInstance,"PSFLSHdll");
		PUREFLDdll = (fp_PUREFLDdllTYPE) GetProcAddress(RefpropdllInstance,"PUREFLDdll");
		QMASSdll = (fp_QMASSdllTYPE) GetProcAddress(RefpropdllInstance,"QMASSdll");
		QMOLEdll = (fp_QMOLEdllTYPE) GetProcAddress(RefpropdllInstance,"QMOLEdll");
		SATDdll = (fp_SATDdllTYPE) GetProcAddress(RefpropdllInstance,"SATDdll");
		SATEdll = (fp_SATEdllTYPE) GetProcAddress(RefpropdllInstance,"SATEdll");
		SATHdll = (fp_SATHdllTYPE) GetProcAddress(RefpropdllInstance,"SATHdll");
		SATPdll = (fp_SATPdllTYPE) GetProcAddress(RefpropdllInstance,"SATPdll");
		SATSdll = (fp_SATSdllTYPE) GetProcAddress(RefpropdllInstance,"SATSdll");
		SATTdll = (fp_SATTdllTYPE) GetProcAddress(RefpropdllInstance,"SATTdll");
		SETAGAdll = (fp_SETAGAdllTYPE) GetProcAddress(RefpropdllInstance,"SETAGAdll");
		SETKTVdll = (fp_SETKTVdllTYPE) GetProcAddress(RefpropdllInstance,"SETKTVdll");
		SETMIXdll = (fp_SETMIXdllTYPE) GetProcAddress(RefpropdllInstance,"SETMIXdll");
		SETMODdll = (fp_SETMODdllTYPE) GetProcAddress(RefpropdllInstance,"SETMODdll");
		SETREFdll = (fp_SETREFdllTYPE) GetProcAddress(RefpropdllInstance,"SETREFdll");
		SETUPdll = (fp_SETUPdllTYPE) GetProcAddress(RefpropdllInstance,"SETUPdll");
		SPECGRdll = (fp_SPECGRdllTYPE) GetProcAddress(RefpropdllInstance,"SPECGRdll");
		SUBLPdll = (fp_SUBLPdllTYPE) GetProcAddress(RefpropdllInstance,"SUBLPdll");
		SUBLTdll = (fp_SUBLTdllTYPE) GetProcAddress(RefpropdllInstance,"SUBLTdll");
		SURFTdll = (fp_SURFTdllTYPE) GetProcAddress(RefpropdllInstance,"SURFTdll");
		SURTENdll = (fp_SURTENdllTYPE) GetProcAddress(RefpropdllInstance,"SURTENdll");
		TDFLSHdll = (fp_TDFLSHdllTYPE) GetProcAddress(RefpropdllInstance,"TDFLSHdll");
		TEFLSHdll = (fp_TEFLSHdllTYPE) GetProcAddress(RefpropdllInstance,"TEFLSHdll");
		THERM0dll = (fp_THERM0dllTYPE) GetProcAddress(RefpropdllInstance,"THERM0dll");
		THERM2dll = (fp_THERM2dllTYPE) GetProcAddress(RefpropdllInstance,"THERM2dll");
		THERM3dll = (fp_THERM3dllTYPE) GetProcAddress(RefpropdllInstance,"THERM3dll");
		THERMdll = (fp_THERMdllTYPE) GetProcAddress(RefpropdllInstance,"THERMdll");
		THFLSHdll = (fp_THFLSHdllTYPE) GetProcAddress(RefpropdllInstance,"THFLSHdll");
		TPFLSHdll = (fp_TPFLSHdllTYPE) GetProcAddress(RefpropdllInstance,"TPFLSHdll");
		TPRHOdll = (fp_TPRHOdllTYPE) GetProcAddress(RefpropdllInstance,"TPRHOdll");
		TQFLSHdll = (fp_TQFLSHdllTYPE) GetProcAddress(RefpropdllInstance,"TQFLSHdll");
		TRNPRPdll = (fp_TRNPRPdllTYPE) GetProcAddress(RefpropdllInstance,"TRNPRPdll");
		TSFLSHdll = (fp_TSFLSHdllTYPE) GetProcAddress(RefpropdllInstance,"TSFLSHdll");
		VIRBdll = (fp_VIRBdllTYPE) GetProcAddress(RefpropdllInstance,"VIRBdll");
		VIRCdll = (fp_VIRCdllTYPE) GetProcAddress(RefpropdllInstance,"VIRCdll");
		WMOLdll = (fp_WMOLdllTYPE) GetProcAddress(RefpropdllInstance,"WMOLdll");
		XMASSdll = (fp_XMASSdllTYPE) GetProcAddress(RefpropdllInstance,"XMASSdll");
		XMOLEdll = (fp_XMOLEdllTYPE) GetProcAddress(RefpropdllInstance,"XMOLEdll");
		
		// If the fluid name starts with the string "REFPROP-", chop off the "REFPROP-"
		if (!strncmp(Ref,"REFPROP-",8))
		{
		char *REFPROPRef=NULL,*RefCopy=NULL;
		double prop;
			
		// Allocate space for refrigerant name
			RefCopy=malloc(strlen(Ref)+1);
		// Make a backup copy
			strcpy(RefCopy,Ref);
		// Chop off the "REFPROP-"
			REFPROPRef = strtok(RefCopy,"-");
			REFPROPRef = strtok(NULL,"-");
		// Run with the stripped Refrigerant name
			prop=REFPROP(Output,Name1,Prop1,Name2,Prop2,REFPROPRef);
		// Free allocated memory
			free(RefCopy);
		// Return the new value
			return prop;
		}
		
		if (!strncmp(Ref,"MIX",3))
		{
			// Sample is "REFPROP-MIX:R32[0.697615]&R125[0.302385]"
			char *REFPROPRef=NULL,*RefCopy=NULL,RefString[255],*Refs[20],*Refrigerant;
			double molefraction;

			//Set global fluid type flag
			FluidType=FLUIDTYPE_REFPROP;
			// Allocate space for refrigerant name
			RefCopy=malloc(strlen(Ref)+1);
			// Make a backup copy
			strcpy(RefCopy,Ref);
			// Chop off the "MIX"
			REFPROPRef = strtok(RefCopy,":");
			i=1;
			while (REFPROPRef!=NULL)
			{
				Refs[i-1]=strtok(NULL,"&");
				if (Refs[i-1]==NULL)
				{
					i--;
					break;
				}
				else
					i++;
			}
			//Flush out RefString
			sprintf(RefString,"");
			for (j=0;j<i;j++)
			{	
				//Get component and its mole fraction
				Refrigerant=strtok(Refs[j],"[]");
				molefraction=strtod(strtok(NULL,"[]"),NULL);
				x[j]=molefraction;
				if (j==0)
					sprintf(RefString,"%s%s.fld",RefString,Refs[j]);
				else
					sprintf(RefString,"%s|%s.fld",RefString,Refs[j]);
			}
			// Free allocated memory
			free(RefCopy);
		}
		else if (!strcmp(Ref,"R508B"))
		{
			i=2;
			strcpy(RefString,"R23.fld|R116.fld");
			x[0]=0.62675;
			x[1]=0.37325;
		}
		else if (!strcmp(Ref,"R410A"))
		{
			i=2;
			strcpy(RefString,"R32.fld|R125.fld");
			x[0]=0.697615;
			x[1]=0.302385;
		}
		else if (!strcmp(Ref,"R404A"))
		{
			i=3;
			strcpy(RefString,"R125.fld|R134a.fld|R143a.fld");
			x[0]=0.35782;
			x[1]=0.038264;
			x[2]=0.60392;
		}
        else if (!strcmp(Ref,"Air"))
		{
			i=3;
			strcpy(RefString,"Nitrogen.fld|Oxygen.fld|Argon.fld");
			x[0]=0.7812;
			x[1]=0.2096;
			x[2]=0.0092;
		}
		else
		{
			i=1;
			strcpy(RefString,"");
			strcat(RefString,Ref);
			strcat(RefString,".fld");
			x[0]=1.0;     //Pure fluid
		}

		strcpy(hf,RefString);
		strcpy(hfmix,"hmx.bnc");
		strcpy(hrf,"DEF");
		strcpy(herr,"Ok");
		
		// If the name of the refrigerant doesn't match 
		// that of the currently loaded refrigerant
		if (strcmp(LoadedREFPROPRef,Ref))
		{
			//...Call SETUP to initialize the program
			SETUPdll(&i, hf, hfmix, hrf, &ierr, herr,
				refpropcharlength*ncmax,refpropcharlength,
				lengthofreference,errormessagelength);
			if (ierr != 0) printf("REFPROP setup gives this error during SETUP: %s\n",herr);
			//Copy the name of the loaded refrigerant back into the temporary holder
			strcpy(LoadedREFPROPRef,Ref);
		}

		// Get the molar mass of the fluid
		WMOLdll(x,&MW);
		if (Output=='B')
		{
			// Critical temperature
			CRITPdll(x,&Tcrit,&pcrit,&dcrit,&ierr,herr,255);
			return Tcrit;
		}
		else if (Output=='E')
		{
			// Critical pressure
			CRITPdll(x,&Tcrit,&pcrit,&dcrit,&ierr,herr,255);
			return pcrit;
		}
		else if (Output=='R')
		{
			long icomp;
			double wmm,Ttriple,tnbpt,tc,pc,Dc,Zc,acf,dip,Rgas;
			// Triple point temperature
			icomp=1;
			if (i>1)
			{
				fprintf(stderr,"Error: Triple point temperature only defined for pure fluids\n");
				return 200;
			}
			INFOdll(&icomp,&wmm,&Ttriple,&tnbpt,&tc,&pc,&Dc,&Zc,&acf,&dip,&Rgas);
			return Ttriple;
		}
		else if (Output=='M')
		{
			// mole mass
			return MW;
		}
		else if (Name1=='T' && Name2=='P')
		{
			// T in K, P in kPa

			// Use flash routine to find properties
			T=Prop1;
			p=Prop2;  
			TPFLSHdll(&T,&p,x,&d,&dl,&dv,xliq,xvap,&q,&e,&h,&s,&cv,&cp,&w,&ierr,herr,errormessagelength);
			if (Output=='H') return h/MW;
			else if (Output=='D') return d*MW;
			else if (Output=='S') return s/MW;
			else if (Output=='U') return e/MW;
			else if (Output=='C') return cp/MW;
			else if (Output=='O') return cv/MW;
			else if (Output=='V') 
			{
				TRNPRPdll(&T,&d,x,&eta,&tcx,&ierr,herr,errormessagelength);
				return eta/1.0e6; //uPa-s to Pa-s
			} 
			else if (Output=='L')
			{
				TRNPRPdll(&T,&d,x,&eta,&tcx,&ierr,herr,errormessagelength);
				return tcx/1000.0; //W/m-K to kW/m-K
			}
			else
				return _HUGE;
		}
		else if (Name1=='T' && Name2=='D')
		{
			// T in K, D in kg/m^3
			// This is the explicit formulation of the EOS
			T=Prop1;
			rho=Prop2/MW;
			
			TDFLSHdll(&T,&rho,x,&p,&dl,&dv,xliq,xvap,&q,&e,&h,&s,&cv,&cp,&w,&ierr,herr,errormessagelength);

			if (Output=='P')
			{
				return p;
			}
			if (Output=='H')
			{
				return h/MW;
			}
			else if (Output=='S')
			{
				return s/MW;
			}
			else if (Output=='U')
			{
				return (h-p/rho)/MW;
			}
			else if (Output=='C')
			{
				return cp/MW;
			}
			else if (Output=='O')
			{
				return cv/MW;
			}
			else if (Output=='V') 
			{
				TRNPRPdll(&T,&rho,x,&eta,&tcx,&ierr,herr,errormessagelength);
				return eta/1.0e6; //uPa-s to Pa-s
			} 
			else if (Output=='L')
			{
				TRNPRPdll(&T,&rho,x,&eta,&tcx,&ierr,herr,errormessagelength);
				return tcx/1000.0; //W/m-K to kW/m-K
			}
			else
				return _HUGE;
		}
		else if (Name1=='T' && Name2=='Q')
		{
			T=Prop1;
			Q=Prop2;
			// Saturation Density
			SATTdll(&T,x,&i,&p,&dl,&dv,xliq,xvap,&ierr,herr,errormessagelength);
			if (Output=='D') 
			{
				return 1/(Q/dv+(1-Q)/dl)*MW;
			}
			else if (Output=='P') 
			{
				PRESSdll(&T,&dl,xliq,&pl);
				PRESSdll(&T,&dv,xvap,&pv);
				return (pv*Q+pl*(1-Q));
			}
			else if (Output=='H') 
			{
				ENTHALdll(&T,&dl,xliq,&hl);
				ENTHALdll(&T,&dv,xvap,&hv);
				return (hv*Q+hl*(1-Q))/MW; // J/kg to kJ/kg
			}
			else if (Output=='S') 
			{
				ENTROdll(&T,&dl,xliq,&sl);
				ENTROdll(&T,&dv,xvap,&sv);
				return (sv*Q+sl*(1-Q))/MW; // J/kg-K to kJ/kg-K
			}
			else if (Output=='U') 
			{
				ENTHALdll(&T,&dl,xliq,&hl);
				ENTHALdll(&T,&dv,xvap,&hv);
				ul=hl-p/dl;
				uv=hv-p/dv;
				return (uv*Q+ul*(1-Q))/MW; // J/kg to kJ/kg
			}
			else if (Output=='C') 
			{
				d=1/(Q/dv+(1-Q)/dl);
				CVCPdll(&T,&d,x,&cv,&cp);
				return cp/MW; // J/kg-K to kJ/kg-K
			}
			else if (Output=='O') 
			{
				d=1/(Q/dv+(1-Q)/dl);
				CVCPdll(&T,&d,x,&cv,&cp);
				return cv/MW; // J/kg-K to kJ/kg-K
			}
			else if (Output=='V') 
			{
				d=1/(Q/dv+(1-Q)/dl);
				TRNPRPdll(&T,&d,x,&eta,&tcx,&ierr,herr,errormessagelength);
				return eta/1.0e6; //uPa-s to Pa-s
			}
			else if (Output=='L') 
			{
				d=1/(Q/dv+(1-Q)/dl);
				TRNPRPdll(&T,&d,x,&eta,&tcx,&ierr,herr,errormessagelength);
				return tcx/1000.0; //W/m-K to kW/m-K
			}
			else
				return _HUGE;
		}
		else if (Name1=='P' && Name2=='Q')
		{
			p=Prop1;
			Q=Prop2;
			// Saturation Density
			SATPdll(&p,x,&i,&T,&dl,&dv,xliq,xvap,&ierr,herr,errormessagelength);
			if (Output=='T')
			{
				return T;
			}
			else if (Output=='D') 
			{
				return 1/(Q/dv+(1-Q)/dl)*MW;
			}
			else if (Output=='P') 
			{
				PRESSdll(&T,&dl,xliq,&pl);
				PRESSdll(&T,&dv,xvap,&pv);
				return (pv*Q+pl*(1-Q));
			}
			else if (Output=='H') 
			{
				ENTHALdll(&T,&dl,xliq,&hl);
				ENTHALdll(&T,&dv,xvap,&hv);
				return (hv*Q+hl*(1-Q))/MW; // J/kg to kJ/kg
			}
			else if (Output=='S') 
			{
				ENTROdll(&T,&dl,xliq,&sl);
				ENTROdll(&T,&dv,xvap,&sv);
				return (sv*Q+sl*(1-Q))/MW; // J/kg-K to kJ/kg-K
			}
			else if (Output=='U') 
			{
				ENTHALdll(&T,&dl,xliq,&hl);
				ENTHALdll(&T,&dv,xvap,&hv);
				ul=hl-p/dl;
				uv=hv-p/dv;
				return (uv*Q+ul*(1-Q))/MW; // J/kg to kJ/kg
			}
			else if (Output=='C') 
			{
				d=1/(Q/dv+(1-Q)/dl);
				CVCPdll(&T,&d,x,&cv,&cp);
				return cp/MW; // J/kg-K to kJ/kg-K
			}
			else if (Output=='O') 
			{
				d=1/(Q/dv+(1-Q)/dl);
				CVCPdll(&T,&d,x,&cv,&cp);
				return cv/MW; // J/kg-K to kJ/kg-K
			}
			else if (Output=='V') 
			{
				d=1/(Q/dv+(1-Q)/dl);
				TRNPRPdll(&T,&d,x,&eta,&tcx,&ierr,herr,errormessagelength);
				return eta/1.0e6; //uPa-s to Pa-s
			}
			else if (Output=='L') 
			{
				d=1/(Q/dv+(1-Q)/dl);
				TRNPRPdll(&T,&d,x,&eta,&tcx,&ierr,herr,errormessagelength);
				return tcx/1000.0; //W/m-K to kW/m-K
			}
			else
				return _HUGE;
		}
		else
			return _HUGE;
}
#endif