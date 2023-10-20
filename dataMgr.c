// ##################################################
// #   dataMgr.c - finally revised on Dec 2018      #
// #   coded by Taeyoon Kim                         #
// #   Copyright (C) 2005 - 2018, Taeyoon Kim,      #
// #   All rights reserved.                         #
// ##################################################
// This file contains functions performing the modification or loading of 
// data.

/*------------------------ Related to network data ---------------------------*/

void ExtractConfigSubroutine(double *r, double *cBnd, double cNeiEdge, 
		int m, int *CS, int *CS2) {
  int k, CS3;

  ConvertRectDomainVector(r, 1);
  FOR_NDIM(k) {
	CS3 = 1;
	if (pbc[k] == 1 || (pbc[k] != 1 && iCell[k] > 0 
			&& iCell[k] < nCell[k] - 1)) {
		if (r[k] >= P2A(bnd.r,0,k,NDIM) && r[k] < P2A(bnd.r,1,k,NDIM)) { 
			(*CS)++; 
			CS3 = 0;
		}
	}
	else if ((iCell[k] == 0 && r[k] < P2A(bnd.r,1,k,NDIM)) 
			|| (iCell[k] == nCell[k] - 1 && r[k] >= P2A(bnd.r,0,k,NDIM))) {
		(*CS)++;
		CS3 = 0;
	}
	if (m == 1 && CS3 == 1 && nCell[k] > 1 
			&& ((r[k] >= P2A(cBnd,0,k,NDIM) - cNeiEdge
			&& r[k] < P2A(cBnd,0,k,NDIM)) || (r[k] >= P2A(cBnd,1,k,NDIM)  
			&& r[k] < P2A(cBnd,1,k,NDIM) + cNeiEdge))) {
		(*CS2)++;  
	}
  }
}

// After the data of a whole network are loaded at LoadConfig(), information
// of particles belonging to the current subdomain is extracted from the data
// here. The rest of data are discarded after this function.
void ExtractConfig(void) {
  int m, n, k, CS, CS2, ind, curr, nIFilaP, nActM, nAbpM, cnt, *pArr;
  double rBnd[NDIM*2], cNeiEdge, len;
  ListInt iFila;


  MALLOC(iFila.l, int, nAct);
  memset(iAct, -1, sizeof(int) * nAct);
  memset(iAbp, -1, sizeof(int) * nAbp);
  cNeiEdge = (motWalk.gTgl != 0) ? neiEdge + 1. : neiEdge;
  V6COPY(rBnd,bnd.r);
  FOR_NDIM(k) {
	CONT(!(pbc[k] == 1));
	if (iCell[k] == 0 && nCell[k] > 1) 
	{ P2A(rBnd,0,k,NDIM) = rGrid[k][nGrid[k] - 1]; }
	else if (iCell[k] == nCell[k] - 1 && nCell[k] > 1) 
	{ P2A(rBnd,1,k,NDIM) = rGrid[k][0]; }
  }

  nActM = 0;
  memset(iFila.l, -1, sizeof(int) * nAct);
  iFila.c = 0;
  FOR_ACT(n) {
	pArr = &P2A(chAct,n,0,nChAc);
    // Count the number of free actins
    if (pArr[0] < 0 && pArr[1] < 0) {
        nActM++;
    }
    // "iFila.l" is the index of actin filaments.
    if (pArr[1] < 0 && pArr[0] > -1) {
        curr = n;
        while(P2A(chAct,curr,0,nChAc) > -1) {
            iFila.l[curr] = iFila.c;
            curr = P2A(chAct,curr,0,nChAc);
        }
        iFila.l[curr] = iFila.c;
        (iFila.c)++;
    }
  }
  if (actSev.gTgl != 0 || actNuc.gTgl != 0 || actAnn.gTgl != 0) {
	nIFilaP = (int)((nActGoal - iFila.c) / nCpu);
	for(n = 0; n < nIFilaP; n++) {
		iFilaP.l[n] = iFila.c + rank * nIFilaP + n;
	}
	iFilaP.c = nIFilaP;
  }
  nAbpM = 0;
  FOR_ABP(n) {
   	pArr = &P2A(chAbp,n,0,nChAb);
    // Count the number of free ABPs
    if (ISMTF(pArr[2])) {
		if (pArr[0] < 0 && pArr[1] < 0 && pArr[3] < 0 && pArr[4] < 0) {
        	nAbpM++;
		}
	}
	else {
		if (pArr[0] < 0 && pArr[1] < 0) {
        	nAbpM++;
		}
	}
  }
  if (nActM % nCpu == 0) { nActM = nActM / nCpu; }
  else { nActM = (int)(nActM / nCpu) + 1; }
  if (nAbpM % nCpu == 0) { nAbpM = nAbpM / nCpu; }
  else { nAbpM = (int)(nAbpM / nCpu) + 1; }

  nActFilaMe = 0;
  nActMe = 0;
  nActCp = 0;
  nAbpMe = 0;
  nAbpCp = 0;
  // If actin is located within the current subdomain, they are listed
  // in iAct and act.id. For example, the actin having the absolute index of 
  // 10 is found for the first time as a particle belonging to the current
  // subdomain, iAct[0] = 10 and act.id[10] = 0. In other words,  
  // act.id[relative index] = absolute index
  // iAct[absolute index] = relative index
  // Also, rAct and fixAct are copied to local variables, act.r and act.fix.
  for(m = 0; m < 2; m++) {
	FOR_ACT(n) {
        CONT(P2A(chAct,n,0,nChAc) < 0 && P2A(chAct,n,1,nChAc) < 0);
		CS = 0; 
		CS2 = 0;
		ExtractConfigSubroutine(&P2(rAct,n,0), rBnd, cNeiEdge, m, &CS, &CS2);
		if ((m == 0 && CS == NDIM) || 
				(m == 1 && CS + CS2 == NDIM && CS != NDIM)) {
			if (m == 0) { ind = nActMe; }
			else { ind = nActMe + nActCp; }
			V3COPY(&P2(act.r,ind,0), &P2(rAct,n,0));
			V3SET_ALL(&P2(act.fBr,ind,0), 0.);
			act.id[ind] = n;
			act.iF[ind] = iFila.l[n];
			iAct[n] = ind;
			if (m == 0) { 
				nActMe++; 
				if (fixAct[n] > -1) { 
					act.fix[ind] = fixAct[n]; 
				}
                if (bndMat.gTgl != 0 && bndUnb.gTgl != 0) {
    	        	act.nFA[ind] = (fixAct[n] > -1) ? 1 : 0;
	            }
			}
			else { nActCp++; }
		}	
	}
    CONT(m == 1);
    cnt = 0;
    FOR_ACT(n) {
		CONT(!(P2A(chAct,n,0,nChAc) < 0 && P2A(chAct,n,1,nChAc) < 0));
		if (cnt >= nActM * rank && cnt < nActM * (rank + 1)) {
			iAct[n] = nActMe;
			UpdateActinAbpMonomerListSubroutine2(iAct[n], n, 0);
			InsertElement1dArrayWoChk(actM.l, &actM.c, n);
			nActMe++;
		}
		cnt++;
    }
  }
  // Perform the same thing as above for ABPs.
  for(m = 0; m < 2; m++) {
	for(n = 0; n < nAbp; n++) {
   		pArr = &P2A(chAbp,n,0,nChAb);
        CONT((ISMTF(pArr[2]) && pArr[0] < 0 && pArr[1] < 0 && pArr[3] < 0 
				&& pArr[4] < 0) || (!(ISMTF(pArr[2])) && pArr[0] < 0 
				&& pArr[1] < 0));
		CS = 0; 
		CS2 = 0;
		ExtractConfigSubroutine(&P2(rAbp,n,0), rBnd, cNeiEdge, m, &CS, &CS2);
		if ((m == 0 && CS == NDIM) || 
					(m == 1 && CS + CS2 == NDIM && CS != NDIM)) {
			if (m == 0) { ind = nAbpMe; }
			else { ind = nAbpMe + nAbpCp; }
			V3COPY(&P2(abp.r,ind,0), &P2(rAbp,n,0));
			V3SET_ALL(&P2(abp.fBr,ind,0), 0.);
			abp.id[ind] = n;

if (tglRecAbpTurn != 0) {
	SetAllValue1dArrayDouble(&P2A(abpTurn,ind,0,7), 7, -1.);
}

			iAbp[n] = ind;
			if (m == 0) { 
				nAbpMe++; 
				if (motSA.gTgl != 0) { abp.mId[ind] = abpMotId[n]; }
			}
			else { nAbpCp++; }
		}
	}
    CONT(m == 1);
    cnt = 0;
    FOR_ABP(n) {
    	pArr = &P2A(chAbp,n,0,nChAb);
        CONT(!((ISMTF(pArr[2]) && pArr[0] < 0 && pArr[1] < 0 && pArr[3] < 0 
				&& pArr[4] < 0) || (!(ISMTF(pArr[2])) && pArr[0] < 0 
				&& pArr[1] < 0)));
		if (cnt >= nAbpM * rank && cnt < nAbpM * (rank + 1)) {
			iAbp[n] = nAbpMe;
			UpdateActinAbpMonomerListSubroutine2(iAbp[n], n, pArr[2] + 1);
			if (pArr[2] == 2 && gTglImpMotM != 0) {
				InsertElement1dArrayWoChk(motM.l, &motM.c, n);
			}
			if (pArr[2] != 2 && gTglImpAcpM != 0) {
				InsertElement1dArrayWoChk(acpM.l, &acpM.c, n);
			}
			if ((pArr[2] == 2 && gTglImpMotM == 0) 
					|| (pArr[2] != 2 && gTglImpAcpM == 0)) {
				if (P2(rAbp,n,0) == 0. && P2(rAbp,n,1) == 0. 
						&& P2(rAbp,n,2) == 0.) {
					GenRandPosSubdom(&P2(abp.r,iAbp[n],0));
				}
				else {
					V3COPY(&P2(abp.r,iAbp[n],0), &P2(rAbp,n,0));
				}
			}
			nAbpMe++;
		}
		cnt++;
    }
  }
  if (nActMe > nActC) {
	Printf0("Error: the value of nActC is too small!!\n\n");
	exit(-1);
  }
  if (nAbpMe > nAbpC) {
	Printf0("Error: the value of nAbpC is too small!!\n\n");
	exit(-1);
  }
  nActMeCp = nActMe + nActCp;
  nAbpMeCp = nAbpMe + nAbpCp;
  // Copy chain information from global variables to local variables
  FOR_ACT(n) {
	pArr = &P2A(chAct,n,0,nChAc);
    CONT(iAct[n] < 0 || (pArr[0] < 0 && pArr[1] < 0));
    if (pArr[0] > -1 && pArr[1] < 0 && iAct[n] < nActMe)
    { nActFilaMe++; }
	Copy1dArrayInt(&P2A(act.ch,iAct[n],0,nChAc), pArr, nChAc);
  }
  FOR_ABP(n) {
    pArr = &P2A(chAbp,n,0,nChAb);
    CONT(iAbp[n] < 0);
	CONT((ISMTF(pArr[2]) && pArr[0] < 0 && pArr[1] < 0 && pArr[3] < 0 
			&& pArr[4] < 0) || (!(ISMTF(pArr[2])) && pArr[0] < 0 
			&& pArr[1] < 0));
	Copy1dArrayInt(&P2A(abp.ch,iAbp[n],0,nChAb), pArr, nChAb);
  }
  if (rheoWay > 0) {
	meaStreParMe.c = 0;
	for(n = 0; n < meaStrePar.c; n++) {
		ind = meaStrePar.l[n];
		if (iAct[ind] > -1 && iAct[ind] < nActMe) {
			meaStreParMe.l[meaStreParMe.c] = ind;
			meaStreParMe.c++;
		}
	}
	appStraParMe.c = 0;
	for(n = 0; n < appStraPar.c; n++) {
		ind = appStraPar.l[n];
		if (iAct[ind] > -1 && iAct[ind] < nActMe) {
			appStraParMe.l[appStraParMe.c] = ind;
			appStraParMe.c++;
		}
	}
  }
  // Free unnecessary global arrays.
  free(rAct);
  free(rAbp);
  free(chAct);
  free(chAbp);
  free(fixAct);
  if (motSA.gTgl != 0) { free(abpMotId); }
  if (recTraj.gTglCho == 0) {
	free(recTraj.act.l); 
	free(recTraj.abp.l); 
  }
  if (rheoWay > 0) {  
	free(appStraPar.l);
	free(meaStrePar.l); 
  }
  free(iFila.l);
  UpdateNeighborList();
  UpdateChainList();

}

// Load information about network configuration from "Config".
void LoadConfig(int mode) {
  int n, k, l, ind, tempInt, tempInt2[2], tglTraj, tglFix;
  char ch[200];
  FILE *fIn; 

  nMot = 0;
  // mode 0 is for loading data to resume simulation in case of errors
  // mode 1 is for the initial loading of configuration file
  if ((fIn = fopen("Config", "r")) == NULL) {
	Printf0("File doesn't exist: Config\n");
	exit(-1);
  }
  // nAct, nAbp, are dimDom[] are loaded at LoadInitParameter(), not here. 
  while (strcmp(ch, "## Position for actin ##\n")) {
	fgets(ch, 200, fIn); 
  }

  // Positions of actin and ABP
  FOR_ACT(n) {
    fscanf(fIn, "%d\t%lf\t%lf\t%lf", &tempInt, &P2(rAct,n,0), &P2(rAct,n,1),
		&P2(rAct,n,2));
  }

  fscanf(fIn, "\n## Position for ABP ##\n");
  FOR_ABP(n) {
	fscanf(fIn, "%d\t%lf\t%lf\t%lf", &tempInt, &P2(rAbp,n,0), &P2(rAbp,n,1), 
		&P2(rAbp,n,2));
  }

  // Chain information of actin and ABP
  memset(chAct, -1, sizeof(int) * nAct * nChAc);
  fscanf(fIn, "\n## Chain for actin ##\n");
  FOR_ACT(n) {
	fscanf(fIn, "%d", &tempInt);
	for(k = 0; k < 2; k++) {
		fscanf(fIn, "%d", &P2A(chAct,n,k,nChAc));
	}
	for(k = 0; k < confNChAcX; k++) {
		for(l = 0; l < confNChAcY; l++) {
			ind = 2 + k * nChAcY * (nChAcX / confNChAcX) + l;
			fscanf(fIn, "%d", &P2A(chAct,n,ind,nChAc));
		}
	}
  }

  fscanf(fIn, "\n## Chain for ABP ##\n");
  FOR_ABP(n) {
	fscanf(fIn, "%d", &tempInt);
	for(k = 0; k < nChAb; k++) {
		fscanf(fIn, "%d", &P2A(chAbp,n,k,nChAb));
	}
	if (P2A(chAbp,n,2,nChAb) == 2) { nMot++; }
  }

  fscanf(fIn, "\n## Position for membrane ##\n");
  fscanf(fIn, "\n## Chain for membrane ##\n");
  // Related to rheological measuremenets
  tglTraj = (recTraj.gTglCho == 0) ? 1 : 0;
  tglFix = (rheoWay > 0 && gTglLoadNetDataFix != 0) ? 1 : 0;

  while(!(ch[0] == 's' && ch[1] == 't')) {
	fgets(ch, 200, fIn);
  }
  fgets(ch, 200, fIn);
  fscanf(fIn, "meaStrePar = %d\n", &meaStrePar.c);
  for(n = 0; n < meaStrePar.c; n++) {
	fscanf(fIn, "%d\n", &tempInt);
	if (tglFix != 0) { meaStrePar.l[n] = tempInt; }
  }
  fscanf(fIn, "\nappStraPar = %d\n", &appStraPar.c);
  for(n = 0; n < appStraPar.c; n++) {
	fscanf(fIn, "%d\n", &tempInt);
	if (tglFix != 0) { appStraPar.l[n] = tempInt; }
  }
  fscanf(fIn, "\nactTraj = %d\n", &recTraj.act.c);
  if (tglTraj != 0) { recTraj.nActL = recTraj.act.c; }
  for (n = 0 ; n < recTraj.act.c ; n++) {
	fscanf(fIn, "%d\n", &tempInt);
	if (tglTraj != 0) { recTraj.act.l[n] = tempInt; }
  }
  fscanf(fIn, "\nabpTraj = %d\n", &recTraj.abp.c);
  if (tglTraj != 0) { recTraj.nAbpL = recTraj.abp.c; }
  for (n = 0 ; n < recTraj.abp.c ; n++) {
	fscanf(fIn, "%d\n", &tempInt);
	if (tglTraj != 0) { recTraj.abp.l[n] = tempInt; }
  }
  fscanf(fIn, "\nfixAct = %d\n", &tempInt);
  for(n = 0; n < tempInt; n++) {
	fscanf(fIn, "%d\t%d\n", &tempInt2[0], &tempInt2[1]);
	if (gTglLoadNetDataFix != 0) {
		fixAct[tempInt2[0]] = tempInt2[1];
		if (fixAct[tempInt2[0]] == 4 && rheoWay > 0) {
		    fixAct[tempInt2[0]] = 14;
		    InsertElement1dArrayWoChk(meaStrePar.l, &meaStrePar.c, tempInt2[0]);
		    InsertElement1dArrayWoChk(appStraPar.l, &appStraPar.c, tempInt2[0]);
		}
	}
  }
  fclose(fIn);
}

// Subroutine for LoadInitParameter(), which judges "yes" or "no" from 
// the file, "condition".
void CheckAnswerYesNo(char *tag, const char *which, int *tgl) {
  if (strcmp(tag, "yes") == 0) { *tgl = 1; }
  else if (strcmp(tag, "no") == 0) { *tgl = 0; } 
  else { 
	Printf0("Error: only yes or no is allowed (%s).\n\n", which); 
	exit(-1); 
  }
}

/*------------------------ Related to network data ---------------------------*/

/*----------------------- Loading initial parameters -------------------------*/

void LoadInitParameterSubroutine2(char *str, int ind, int *tgl1, int *tgl2) {
  *tgl1 = (ind > 0) ? 1 : 0;
  *tgl2 = (ind == 2) ? 1 : 0;
  if (ind < 0 || ind > 2) {
	Printf0("Error: choice of %s should be either of 0, 1, or 2.\n\n", str);
  	MPI_Barrier(MPI_COMM_WORLD);
	exit(-1);
  }
}

void LoadInitParameterSubroutine(FILE *fIn, const char *which, int *tgl) {
  char tag[80];

  fscanf(fIn, "%s", tag);
  if (tag[0] == 'y' && tag[1] == 'e' && tag[2] == 's') 
  { *tgl = 1; }
  else if (tag[0] == 'n' && tag[1] == 'o') 
  { *tgl = 0; }
  else {
	Printf0("Error: only yes or no is allowed (%s).\n\n", which); 
	exit(-1); 
  }
}

// Load parameters from "Config" and "condition" at the very beginning.
// Those parameters are needed to define arrays, and so on.
void LoadInitParameter(void) {
  int n, k, tempInt[4], CS;
  double tempDbl[3];
  char tag[200], direc[4]= "xyz";
  FILE *fIn;

  nChAb = 5;
  MALLOC(bndUnb.facK0,double,NDIM*2);
  MALLOC(bndUnb.facX,double,NDIM*2);
  MALLOC(bndReb.facK,double,NDIM*2);
  MALLOC(bnd.drag,double,NDIM*2);
  MALLOC(bnd.stfSpr,double,NDIM*2);
  MALLOC(abpF.a90,Force,NK_ABP);
  MALLOC(abpF.cr,Force,NK_ABP);
  MALLOC(abpF.bend,Force,NK_ABP);
  MALLOC(abpF.spr,Force,NK_ABP);

  // From "condition", several types of informatio are loaded: MPI method,
  // rheological methods, dynamic behaviors of ABPs
  if ((fIn = fopen("condition", "r")) == NULL) {
	Printf0("File doesn't exist: condition\n");
	exit(-1); 
  }
  fgets(tag, 200, fIn);

  mpiMethod = 0;
  updSubdSize.prdR = 1.;

  nActPerSeg = 71;
  nChAcX = 10;
  nChAcY = 2;

  fscanf(fIn, "Simulation mode(0: matrix assembly, 1: measurement) = %d\n", 
		&rheoWay);
  if (rheoWay == 0) { rheoWay = -1; }
  else if (rheoWay == 1) { rheoWay = 1; }
 
  if (rheoWay > 0) {
	 gTglLoadNetData = 1;
  }

  dir2D = -1;
  gTglLoadNetDataFix = 1;

  fscanf(fIn, "Duration of matrix formation(s) = %lf\n", &netForm.durR);
  fscanf(fIn, "Duration of measurement(s) = %lf\n", &rheo.durR);


  motActiv.durR = 0.;
  // dirStr: direction of stress/strain
  // dirNoPBC: direction of no periodic boundary condition for bulk rheology
  // dirOther: the other direction
  dirStr = -1;
  dirNoPBC = -1;
  dirOther = -1;


  gTglBead = 0;
  bead.n = 0;
  bead.rad2 = 0.5;
  bead.facStfRep = 1.;
  bead.thkRep = 0.05;
  recBeadLoc.prdR = 0.01;
  beadBind.gTgl = 0;
  bead.rheoWay = 1;

  bead.rheoDir = 0;
  bead.rheoSig = 0;
  bead.rheoPreMag = 0;
  bead.rheoPreRate = 10;
  bead.rheoAmp = 10;
  bead.rheoPrdR = 1.;

  if (rheoWay < 0) {
	gTglActTherm = 0;
	actNuc.gTgl = 1;
	actBch.gTgl = 1;
	actAss.gTgl = 1;
	gTglAcpTherm = 0;
  }
  else {
	gTglActTherm = 1;
	actNuc.gTgl = 0;
	actBch.gTgl = 0;
	actAss.gTgl = 0;
	gTglAcpTherm = 1;
  }
  actNuc.gTglFN = 0;
  gTglActDynNF = 1;
  gTglActDynPres = 1;
  actDis.gTgl = 0;
  actBst.gTgl = 0;
  actCap.gTgl = 0;
  actUnc.gTgl = 0;
  actSev.gTgl = 0;
  actSev.gTglCap[0] = 0;
  actSev.gTglCap[1] = 0;
  actSev.gTglBst[0] = 0;
  actSev.gTglBst[1] = 0;
  actAnn.gTgl = 0;
  actDgd.gTgl = 0;
  gTglAcpDynNF = 1;
  gTglAcpDynPres = 1;
  acpInaUnb.gTgl = 1;
  acpMoBind.gTgl = 1;
  acpUnb.gTgl = 1;
  acpReb.gTgl = 1;
  acpReb.gTglCrsAng = 0;
  gTglImpAcpM = 1;

  gTglMotTherm = 1;
  gTglMotUnbRebNF = 1;
  gTglMotUnbRebPres = 1;
  gTglMotWalkNF = 0;
  gTglMotWalkPres = 1;
  motWalk.gTgl = 1;
  gTglMotWalkSld = 0;
  motSA.gTgl = 0;
  motSA.gTglConSiz = 0;
  motInaUnb.gTgl = 1;
  motMoBind.gTgl = 1;
  motUnb.gTgl = 1;
  motReb.gTgl = 1;
  motReb.gTglCrsAng = 0;
  gTglImpMotM = 0;
  motInaUnb.gTgl = 1;
  motMoBind.gTgl = 1;
  motReb.gTglCrsAng = 1;
  motReb.gTglOppDir = 0;
  motSA.to.gTgl = 0;

  updMono.prdR = 0.1;
  // Parameters for dynamic behaviors of actin and ABP

  actBch.k = 100.;
  actNuc.k = 0.00001;
  actAss.k[0] = 6000;
  actAss.k[1] = 0;
  actDis.k[0] = 0;
  actDis.k[1] = 0;
  actDis.facKWA = 1.0;
  actBst.k[0] = 0.1;
  actBst.k[1] = 0.1;
  actBst.facKWA = 1.0;
  actSev.k = 0;
  actSev.facX = 1.6;
  actSev.facKWA = 1.0;
  actAnn.k = 1.0;
  actAnn.ang = 10.0;
  actCap.k[0] = 1.0;
  actCap.k[1] = 1.0;
  actUnc.k[0] = 1.0;
  actUnc.k[1] = 1.0;
  actDgd.k = 10;
  actDgd.x1 = 0;
  actDgd.x2 = 0;
  actDgd.dist = 150;

  acpReb.facK = 1.0;
  acpMoBind.facK = acpReb.facK;
  acpUnb.facK0 = 0.01;
  acpUnb.facX = 1.0e-10;
  acpUnb.facK0c = 0.1;
  acpUnb.facXc = -1.0e-9;

  motReb.facK = 1.0;
  motUnb.facK0 = 1.0;
  motMC.nHead = 4;
  motMC.k01 = 40;
  motMC.k10 = 2;
  motMC.k12 = 1000;
  motMC.k21 = 1000;
  motMC.k20 = 20;
  motSA.nMotPerTF = 1;
  updMotN.prdR = -1;
  motSA.kAss = 1e-2;
  motSA.to.k = 1e-3;
  motSA.to.k0 = 1e-3;
  motSA.to.x = -1e-10;

  netForm.facK = 100;
  pres.facK = 1;
  motActiv.facK = 1;
  actF.rep.facStf = 1.0;
  actF.bend.facStf = 2.24;
  actF.spr.facStf = 10;
  abpF.rep.facStf = 0;
  abpF.bend[0].facStf = 1.0;
  abpF.bend[1].facStf = 1.0;
  abpF.bend[2].facStf = 1.0;
  abpF.a90[0].facStf = 1.0;
  abpF.a90[1].facStf = 0;
  abpF.a90[2].facStf = 0;
  abpF.spr[0].facStf = 1.0;
  abpF.spr[1].facStf = 1.0;
  abpF.spr[2].facStf = 1.0;

  bnd.radRnd = -1;
  bndMv.gTgl = 0;
  bndMv.stf[0] = 10000;
  bndMv.stf[3] = 10000;
  bndMv.stf[1] = 10000;
  bndMv.stf[4] = 10000;
  bndMv.stf[2] = 10000;
  bndMv.stf[5] = 10000;
  bndMv.thk[0] = 5;
  bndMv.thk[3] = 5;
  bndMv.thk[1] = 5;
  bndMv.thk[4] = 5;
  bndMv.thk[2] = 5;
  bndMv.thk[5] = 5;
  bndMv.stfUnit = 0;
  bndVol.gTgl = 0;
  bndVol.facStf = 1;
  bndUnb.gTgl = 0;
  bndUnb.facK0[0] = 1.0;
  bndUnb.facK0[1] = 1.0;
  bndUnb.facK0[2] = 1.0;
  bndUnb.facK0[3] = 1.0;
  bndUnb.facK0[4] = 1.0;
  bndUnb.facK0[5] = 1.0;
  bndUnb.facX[0] = 1.0;
  bndUnb.facX[1] = 1.0;
  bndUnb.facX[2] = 1.0;
  bndUnb.facX[3] = 1.0;
  bndUnb.facX[4] = 1.0;
  bndUnb.facX[5] = 1.0;
  tempInt[0] = 2;
  bndReb.gTglPa = 2;
  bndReb.depR = 200;
  bndReb.facK[0] = 0;
  bndReb.facK[1] = 0;
  bndReb.facK[2] = 1000;
  bndReb.facK[3] = 1000;
  bndReb.facK[4] = 0;
  bndReb.facK[5] = 0;
  bndReb.gTglSpr = 0;
  bnd.stfSpr[0] = 1.0;
  bnd.stfSpr[1] = 1.0;
  bnd.stfSpr[2] = 1.0;
  bnd.stfSpr[3] = 1.0;
  bnd.stfSpr[4] = 1.0;
  bnd.stfSpr[5] = 1.0;
  bnd.gTglActMv = 0;
  bnd.drag[0] = 100;
  bnd.drag[1] = 100;
  bnd.drag[2] = 100;
  bnd.drag[3] = 100;
  bnd.drag[4] = 100;
  bnd.drag[5] = 100;
  tempInt[0] = 0;
  bndMat.gTglPa = 0;
  bndMat.facK0 = 1;
  bndMat.facX = 1;
  bndMat.maxNFA = 20;
  recBndLoc.prdR = -1;
  recBndUnbReb.prdR = 1;
  recBndActMat.prdR = -1;
  recBndTracF.prdR = -1;
  
  gTglMb = 0;
  gTglNuc = 0;
  gTglLoadMbNucData = 0;
  dimMbNuc = 3;
  nObjMbNuc = 1;
  radMb = 3;
  memb.len = 0.35;
  memb.thk = 0.1;
  memb.dragR = 1;
  memb.bend.facStf = 1;
  memb.spr.facStf = 1;
  memb.spr2.facStf = 1;
  mbVol.gTgl = 1;
  memb.vol.facStf = 1;
  mbAre.gTgl = 1;
  memb.area.facStf = 1;

  gTglMbTherm = 1;
  mbPro.gTgl = 0;
  mbPro.facK = 1;
  mbPro.gTglCT = 0;
  mbPro.rCT[0] = 0.5;
  mbPro.rCT[1] = 0.5;
  mbPro.rCT[2] = 0.5;
  mbPro.f = 100;
  mbPro.durR = 5;
  mbFix.gTgl = 0;
  mbFix.facK = 0.01;
  mbDef.gTgl = 0;
  mbDef.facK0 = 1;
  mbDef.facX = 1;
  sideMb = 0;
  thkMbActNuc = 0.1;
  gTglMbCont = 0;
  tempInt[0] = 2;
  mbReb.gTglPa = 0;

  mbReb.por = 1.0;
  nMbAct = 5;
  mbReb.facK = 1;
  mbUnb.gTgl = 0;
  mbUnb.facK0 = 1;
  mbUnb.facX = 1;
  tempInt[0] = 0;
  mbMat.gTglPa = 0;
  mbMat.facK0 = 1;
  mbMat.facX = 8;
  mbMat.maxNFA = 20;

  tempInt[0] = 2;
  mbSld.act.por = 0.1;
  mbSld.abp.por = 0.1;
  mbSld.critDist = 300;
  mbSld.eqDist = 100;
  mbSld.stf = 0.0003;
  mbSld.drag = 10;

  radNuc = 2;
  memb.nucLen = 0.35;
  memb.nucThk = 0.2;
  memb.nucDragR = 1;
  memb.nucBend.facStf = 1;
  memb.nucSpr.facStf = 1;
  memb.nucSpr2.facStf = 1;
  nucVol.gTgl = 1;
  memb.nucVol.facStf = 1;
  nucAre.gTgl = 1;
  memb.nucArea.facStf = 1;
  gTglNucTherm = 0;
  tempInt[0] = 0;
  nucReb.gTglPa = 0;
  nucReb.por = 1;
  nNucAct = 5;
  nucReb.facK = 1;
  nucUnb.gTgl = 0;
  nucUnb.facK0 = 1;
  nucUnb.facX = 1;

  recConf.prdR = 1;
  recConfMlb.prdR = -1;
  recConfVmd.prdR = 1;
  recConfVmd.mode = 1;
  recConfVmd.gTglBnd = 1;
  recConfVmd.gTglInfo = 1;
  recConfVmd.minF = 0;
  recConfVmd.maxF = 100;
  recFilaL.prdR = 1;
  recMotSize.prdR = -1;
  recMotPos.prdR = -1;
  recCrsDist.prdR = 1;
  recPoreSize.prdR = -1;
  recConn.prdR = -1;
  recPerc.prdR = -1 ;
  findSupp.prdR = -1;
  porFindSupp = 0.2;
  kindFindSupp = 0;

  recLongF.prdR = 1;
  recE.prdR = 1;
  recSecStre.prdR = 1;
  nSecStreDiv[0] = 20;
  nSecStreDiv[1] = 20;
  nSecStreDiv[2] = 20;
  recSecStre.mode = 0;

  recActSev.prdR = -1;
  recAbpUnb.prdR = -1;
  recAbpBind.prdR = -1;
  recAbpTurn.prdR = -1;
  recAbpDyn.prdR = -1;

  recMbCen.prdR = -1;
  recMbDim.prdR = -1;
  recInfo.prdR = 1;

  fscanf(fIn, "Period of recording Output and Progress(s) = %lf\n", 
		&recProg.prdR); 
   
  fgets(tag, 200, fIn);
  fscanf(fIn, "Width of domain(x, y, z in um) = %lf, %lf, %lf\n", 
		&tempDbl[0], &tempDbl[1], &tempDbl[2]);
  if (gTglLoadNetData == 0) {
	FOR_NDIM(n) { dimDom[n] = L_UM2S(tempDbl[n]); }
  }
  V3SET(pbc, 1, 0, 0);

  fscanf(fIn, "Fiber concentration(in mg/ml or given) = %s\n", tag);
  if (strcmp(tag, "given") == 0) { 
	if (gTglLoadNetData == 0) {
		Printf0("Error: no given actin concentration without "
				"loaded network data.\n\n"); 
		exit(-1); 
	}
	cAct = -1.; 
  }
  else { cAct = atof(tag); }
  // R values of ACPC, ACPB, and motor
  fscanf(fIn, "Bundler density(R value or given) = %s\n", tag);
  RAbp[0] = (rheoWay > 0) ? -1. : atof(tag);

  fscanf(fIn, "Cross-linker density(R value or given) = %s\n", tag);
  RAbp[1] = (rheoWay > 0) ? -1. : atof(tag);

  RAbp[2] = 0.;

  if (gTglLoadNetData == 0 && (RAbp[0] == -1. || RAbp[1] == -1. 
		|| RAbp[2] == -1.)) {
	Printf0("Error: no given ABP concentration without "
			"loaded network data.\n\n"); 
	exit(-1); 
  }

  fscanf(fIn, "Maximum length of filaments(in um) = %lf\n", &maxFilaLen);
  fscanf(fIn, "Angle between fibers(in deg) = %lf\n", &angFila);
  fscanf(fIn, "Limit of relative position of branching points(0-1) = %lf, "
		"%lf\n", &limPosBch[0], &limPosBch[1]);

  fgets(tag, 200, fIn);
  fscanf(fIn, "Magnitude of strain(unitless) = %lf\n", &pres.mag);
  fscanf(fIn, "Strain rate(1/s) = %lf\n", &pres.rate);
  fscanf(fIn, "Period of recording stress and strain(s) = %lf\n", 
		&recStre.prdR);
  bulkRheoWay = 1;
  bulkRheoType = 0;
  signStr = 1;
  dirStr = 0;
  dirNoPBC = 1;
  sinuStr.amp = 0.0000001;
  sinuStr.prdR = 1000.;
  stra.lim[0] = -10;
  stra.lim[1] = 3;

  fclose(fIn);

  if (cAct >= 0. && actAss.gTgl == 0 && actDis.gTgl == 0 
		&& actBst.gTgl == 0) { 
	cAct = -1.; 
  }
  if (gTglLoadNetData == 0 && actAss.gTgl == 0 && actNuc.gTgl == 0) {
	Printf0("Error: actin assembly and nucleation should be allowed "
			"at least.\n\n"); 
	exit(-1); 
  }
  if (rheoWay > 0 && bndMv.gTgl != 0) {
	Printf0("Warning: bulk rheology cannot be used with moving boundaries. "
			"The moving boundaries will be deactivated!!\n\n");
	bndMv.gTgl = 0;
  }
  if (rheoWay == 0 || rheoWay == 2) {
	if ((bndUnb.gTgl != 0 || bndReb.gTgl != 0 || bnd.gTglActMv != 0) 
			&& pbc[0] != 0 && pbc[1] != 0 && pbc[2] != 0) {
		Printf0("Warning: periodic boundary conditions exist in all directions,"
				" so the unbinding/binding of actins on boundaries will be"
				" ignored!\n\n");
		bndUnb.gTgl = bndReb.gTgl = bnd.gTglActMv = 0;
	}
  }

  V3SET_ALL(nAbpDet, 0);
  if (gTglLoadNetData != 0) {
	// From "Config", the values of nAct, nAbp, and dimDom[] are loaded.
	if ((fIn = fopen("Config", "r")) == NULL) {
	    Printf0("File doesn't exist: Config\n");
	    exit(-1);
	}
	fscanf(fIn, "nAct : %d\n", &nAct);
	fscanf(fIn, "nAbp : %d\n", &nAbp);
	fscanf(fIn, "nMb : %d\n", &nMb);
	fscanf(fIn, "dimDom : %lf, %lf, %lf\n", &dimDom[0], &dimDom[1], &dimDom[2]);

	fscanf(fIn, "nActPerSeg, nChAcX, nChAcY : %d, %d, %d\n", &tempInt[0], 
			&confNChAcX, &confNChAcY);
	if (tempInt[0] != nActPerSeg){
	    Printf0("Error: the length of actin cylindrical segments "
				"doesn't match between 'condition' and 'Config': "
				"%g nm in 'condition' vs %g nm in 'Config'\n", 
				nActPerSeg * 7., tempInt[0] * 7.);
	    exit(-1);
	}
	if (confNChAcX > nChAcX) {
		Printf0("Error: the number of binding sites in a longitudinal direction"
				" on each actin segment in 'Config' is greater than that in"
				" 'condition': %d in 'Config' vs %d in 'condition'\n", 
				confNChAcX, nChAcX);
	    exit(-1);
	}
	else if (confNChAcX < nChAcX) {
		if (nChAcX % confNChAcX == 0) {
			Printf0("Warning: the number of binding sites in a longitudinal"
					" direction on each actin segment in 'Config' is smaller"
					" than that in 'condition': %d in 'Config' vs %d in"
					" 'condition'\n\n", confNChAcX, nChAcX);
		}
		else {
			Printf0("Error: the number of binding sites in a longitudinal"
					" direction on each actin segment in 'Config' is smaller"
					" than that in 'condition'. In this case, the latter"
					" should be a multiple of the former: %d in 'Config'"
					" vs %d in 'condition'\n", confNChAcX, nChAcX);
		    exit(-1);
		}
	}
	if (confNChAcY > nChAcY) {
		Printf0("Error: the number of binding sites in a transverse direction"
				" on each actin segment in 'Config' is greater than that in"
				" 'condition': %d in 'Config' vs %d in 'condition'\n", 
				confNChAcY, nChAcY);
	    exit(-1);
	}
	else if (confNChAcY < nChAcY) {
		Printf0("Warning: the number of binding sites in a transverse"
				" direction on each actin segment in 'Config' is smaller"
				" than that in 'condition': %d in 'Config' vs %d in"
				" 'condition'\n\n", confNChAcY, nChAcY);
	}

	while (strcmp(tag, "## Chain for actin ##\n")) {
		fgets(tag, 200, fIn); 
	}
	FOR_ACT(n) {
		for(k = 0; k < 3 + confNChAcX * confNChAcY; k++) { 
			fscanf(fIn, "%d", &tempInt[0]);
		}
	}
	while (strcmp(tag, "## Chain for ABP ##\n")) {
		fgets(tag, 200, fIn); 
	}
	FOR_ABP(n) {
		for(k = 0; k < nChAb + 1; k++) { 
			fscanf(fIn, "%d", &tempInt[0]);
			if (k == 3) { nAbpDet[tempInt[0]]++; }
		}
	}
	while(!(tag[0] == 'r' && tag[1] == 'G')) {
		fgets(tag, 200, fIn); 
	}
	fscanf(fIn, "pbc = %d, %d, %d\n", &confPbc[0], &confPbc[1], &confPbc[2]);
  }
  else {
	nAct = 0;
	nAbp = 0;
  }
  CS = 1;
  if (bndUnb.gTgl == 0) {
	free(bndUnb.facK0);
	free(bndUnb.facX);
  }
  if (bndReb.gTgl == 0) {
	free(bndReb.facK);
  }
  if (bnd.gTglActMv == 0) {
	free(bnd.drag);
  }
  gTglActCapAll = (actCap.gTgl != 0 || actUnc.gTgl != 0
		|| (actSev.gTgl != 0 && gTglActSevCap != 0)) ? 1 : 0;
}

/*----------------------- Loading initial parameters -------------------------*/

/*----------------------- Adding and deleting elements -----------------------*/

void AddFreeActinAbpSubroutine(int nMall, int *nAddMe, int *nAddMeSum) {
  int k, mSum, *mAll;
  double volSum, *volAll;

  volMe = HowManyInOutBound();
  MALLOC(volAll, double, nCpu);
  MALLOC(mAll, int, nCpu);
  MPI_Gather(&volMe, 1, MPI_DOUBLE, volAll, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (rank == 0) {
	volSum = SumArrDbl(volAll, nCpu);
	for(k = 0; k < nCpu; k++) {
		mAll[k] = (int)(volAll[k] / volSum * nMall);
	}
	mSum = SumArrInt(mAll, nCpu);
	mAll[nCpu - 1] += nMall - mSum;
  }
  MPI_Bcast(mAll, nCpu, MPI_INT, 0, MPI_COMM_WORLD);
  *nAddMe = mAll[rank];
  *nAddMeSum = SumArrInt(mAll, rank);

  free(volAll);
  free(mAll);
}

void AddFreeAbp(void) {
  int n, kind, CS, sum, nAbpAddSumMe, nAbpMallSum, tempInt, tempInt2;
  int nAbpAddMe[NK_ABP], sum2[NK_ABP], *nAbpAddSumAll;
  int diff, ind, *nNucAll;

  Printf0("%d ACPC, %d ACPB, %d motor are added.\n",
		nAbpMall[0], nAbpMall[1], nAbpMall[2]);
  nMot += nAbpMall[2];
  nAbpMallSum = V3SUM(nAbpMall);

  if (bnd.gTglRnd != 0) {
	if (motSA.gTgl != 0) {
		AddFreeActinAbpSubroutine(motSA.nNucMe, &tempInt, &tempInt2);
		motSA.nNucMe = tempInt;
	}
	for(n = 0; n < NK_ABP; n++) {
		AddFreeActinAbpSubroutine(nAbpMall[n], &nAbpAddMe[n], &sum2[n]);
	}
	sum = SumArrInt(sum2, NK_ABP);
	nAbpAddSumMe = SumArrInt(nAbpAddMe, NK_ABP);
  }
  else {
	for(n = 0; n < NK_ABP; n++) {
		nAbpAddMe[n] = nAbpMall[n] / nCpu
				+ (rank < (nAbpMall[n] % nCpu) ? 1 : 0);
	}
	if (motSA.gTgl != 0) {
		MALLOC(nNucAll, int, nCpu);
		if (rank == 0) {
			SetAllValue1dArrayInt(nNucAll, nCpu, motSA.nNucMe / nCpu);
			for(n = 0; n < motSA.nNucMe % nCpu; n++) {
				while(1) {
					ind = GenRandIntIndex(nCpu);
					BREAK(nNucAll[ind] == motSA.nNucMe / nCpu);
				}
				nNucAll[ind]++;
			}
			if (SumArrInt(nNucAll, nCpu) == 0 && nAbpMall[2] > 0) {
				ind = GenRandIntIndex(nCpu);
				nNucAll[ind] = 1;
			}
		}
		MPI_Bcast(nNucAll, nCpu, MPI_INT, 0, MPI_COMM_WORLD);
		motSA.nNucMe = nNucAll[rank];
		sum = SumArrInt(nNucAll, rank);
		diff = nAbpMall[2] - sum * motSA.nMotPerSide * 2;
		if (diff <= 0) {
			nAbpAddMe[2] = 0;
		}
		else {
			nAbpAddMe[2] = motSA.nNucMe * motSA.nMotPerSide * 2;
			if (diff - motSA.nNucMe * motSA.nMotPerSide * 2
					< motSA.nMotPerSide * 2 && motSA.nNucMe > 0) {
				nAbpAddMe[2] += diff - motSA.nNucMe * motSA.nMotPerSide * 2;
			}
		}
		free(nNucAll);
	}
	nAbpAddSumMe = SumArrInt(nAbpAddMe, NK_ABP);
	MALLOC(nAbpAddSumAll, int, nCpu);
	MPI_Allgather(&nAbpAddSumMe, 1, MPI_INT, nAbpAddSumAll, 1, MPI_INT,
			MPI_COMM_WORLD);
	sum = SumArrInt(nAbpAddSumAll, rank);
	free(nAbpAddSumAll);
  }
  for(n = 0; n < NK_ABP; n++) {
    nAbpDet[n] += nAbpMall[n];
    nAbpMall[n] = 0;
  }
  nAcpMme += nAbpAddMe[0] + nAbpAddMe[1];
  nMotMme += nAbpAddMe[2];
  // Shift the copied ABPs
  for(n = nAbpCp - 1; n >= 0; n--) {
	UpdateActinAbpMonomerListSubroutine(nAbpMe + nAbpAddSumMe + n,
			nAbpMe + n, 1);
  }
  // Determine the kind of free ABP to add
  for(n = 0; n < nAbpAddSumMe; n++) {
	CS = 0;
	while(CS == 0) {
		kind = GenRandIntIndex(NK_ABP);
		if (nAbpAddMe[kind] > 0) {
			nAbpAddMe[kind]--;
			CS = 1;
		}
	}
	// Add the free ABP
	UpdateActinAbpMonomerListSubroutine2(nAbpMe + n, nAbp + sum + n, kind + 1);
	if ((kind != 2 && gTglImpAcpM == 0) || (kind == 2 && gTglImpMotM == 0)) {
		GenRandPosSubdom(&P2(abp.r,nAbpMe + n,0));
	}
	if (kind != 2 && gTglImpAcpM != 0) {
		InsertElement1dArrayWoChk(acpM.l, &acpM.c, nAbp + sum + n);
	}
	if (kind == 2 && gTglImpMotM != 0) {
		InsertElement1dArrayWoChk(motM.l, &motM.c, nAbp + sum + n);
	}
  }
  // All other free ABPs belonging to other submains should be not here
  for (n = 0; n < nAbpMallSum; n++) {
	if (!(n >= sum && n < sum + nAbpAddSumMe)) {
		iAbp[nAbp + n] = -1;
	}
  }
  nAbpMe += nAbpAddSumMe;
  nAbp += nAbpMallSum;
}

void AddFreeActin(void) {
  int n, k, ind[NDIM], nActAddMe, cnt, *cntAll, sum;
  int begin[NDIM], end[NDIM];
  double dr[NDIM], dist, thres[2], r[NDIM];

  if (nActMall > 0) {
	if (bnd.gTglRnd != 0) {
		AddFreeActinAbpSubroutine(nActMall, &nActAddMe, &sum);
	}
	else {
		nActAddMe = nActMall / nCpu;
		sum = rank * nActAddMe;
	}
	for (n = nActCp - 1; n >= 0; n--) {
		UpdateActinAbpMonomerListSubroutine(nActMe + nActAddMe + n,
				nActMe + n, 0);
	}
	for (n = 0; n < nActAddMe; n++) {
		UpdateActinAbpMonomerListSubroutine2(nActMe + n,
				nAct + sum + n, 0);
		actM.l[n + actM.c] = nAct + sum + n;
	}
	for (n = 0; n < nActMall; n++) {
		CONT(n >= sum && n < sum + nActAddMe);
		iAct[nAct + n] = -1;
	}
	actM.c += nActAddMe;
	nActMe += nActAddMe;
	nAct += nActMall;
	if (rank == 0) {
  		Printf0("%d actins are added.\n", nActMall);
	} 
  } 
}

void DeleteAbp(void) {
  int n, k, CS, ind, ind2, side, cnt, nDelSum, *pArr;
  int *abpIndL, *abpKindL, *chAbp2, *abpMotId2, nDel[NK_ABP], nDel2[NK_ABP];
  int cntAbp, cntInactAbp, CS2;
  double *rAbp2;
  ListInt del;

  for(n = 0; n < NK_ABP; n++) {
	if (nAbpMall[n] < 0) {
		nDel[n] = -1 * nAbpMall[n];
		nAbpDet[n] += nAbpMall[n];
		nAbpMall[n] = 0;
	}
	else { nDel[n] = 0; }
  }
  nMot -= nDel[2];

  MALLOC(del.l,int,nAbp);
  MALLOC(chAbp2,int,nAbp*nChAb);
  MALLOC(rAbp2,double,nAbp*NDIM);
  MALLOC(abpIndL,int,nAbp);
  if (motSA.gTgl != 0) {
	MALLOC(abpMotId2,int,nAbp);
  }

  Printf0("%d ACPC, %d ACPB, %d motor are deleted.\n", 
			nDel[0], nDel[1], nDel[2]); 
  nDelSum = V3SUM(nDel);
  del.c = 0;
  if (rank == 0) {

	cntInactAbp = 0;
	FOR_ABP(n) {
	    CONT(P2A(chAbp,n,2,nChAb) != 1);
	    CONT(P2A(chAbp,n,0,nChAb) > -1 && P2A(chAbp,n,1,nChAb) > -1);
	    cntInactAbp++;
	}
	CS2 = 0;
	cntAbp = 0;

    MALLOC(abpKindL,int,nAbp);
	FOR_ABP(n) { abpKindL[n] = P2A(chAbp,n,2,nChAb); }
	for(n = 0; n < NK_ABP; n++) {
		for(k = 0; k < nDel[n]; k++) {
			CS = 0;
			while (CS == 0) {
				ind = GenRandIntIndex(nAbp);
				CONT(!(abpKindL[ind] == n));
				if (ISMTF(n)) {
					CONT(P2A(chAbp,ind,3,nChAb) > -1 
							&& P2A(chAbp,ind,4,nChAb) > -1);
				}

				CONT(CS2 == 0 && P2A(chAbp,ind,0,nChAb) > -1 
						&& P2A(chAbp,ind,1,nChAb) > -1);
				if (n == 1) {
				    cntAbp++;
				}
				if (cntAbp >= cntInactAbp) { CS2 = 1; }

				del.l[del.c] = ind;
				(del.c)++;
				abpKindL[ind] = -1;
				CS = 1;
			}
		}
	}
	free(abpKindL);
  }
  MPI_Bcast(del.l, nDelSum, MPI_INT, 0, MPI_COMM_WORLD);

  for(n = 0; n < nDelSum; n++) {
	ind = del.l[n];
	for(k = 0; k < 5; k++) {
		CONT(k == 2);
		ind2 = P2A(chAbp,ind,k,nChAb);
		CONT(!(ind2 > -1));
		if (k < 2) {
			side = FindAbpActinChain(ind2, ind, 1);
			P2A(chAct,ind2,side,nChAc) = -1;	
		}
		else {
			P2A(chAbp,ind2,7 - k,nChAb) = -1;	
		}
	}
	SetAllValue1dArrayInt(&P2A(chAbp,ind,0,nChAb), nChAb, -1);
  }

  cnt = 0;
  FOR_ABP(n) {
	CONT(!(P2A(chAbp,n,2,nChAb) > -1));
	Copy1dArrayInt(&P2A(chAbp2,cnt,0,nChAb), &P2A(chAbp,n,0,nChAb), nChAb);
	V3COPY(&P2(rAbp2,cnt,0), &P2(rAbp,n,0));
	if (motSA.gTgl != 0) {
		abpMotId2[cnt] = abpMotId[n];
	}
	for(k = 0; k < 2; k++) {
		ind = P2A(chAbp,n,k,nChAb);
		CONT(!(ind > -1));
		side = FindAbpActinChain(ind, n, 1);
		P2A(chAct,ind,side,nChAc) = cnt;
	}
	abpIndL[n] = cnt;
	cnt++;
  }
  nAbp = cnt;
  FOR_ABP(n) {
	for(k = 3; k < 5; k++) {
		CONT(P2A(chAbp2,n,k,nChAb) < 0);
		P2A(chAbp2,n,k,nChAb) = abpIndL[P2A(chAbp2,n,k,nChAb)];
	}
  }

  FOR_ABP(n) {
    Copy1dArrayInt(&P2A(chAbp,n,0,nChAb), &P2A(chAbp2,n,0,nChAb), nChAb);
    V3COPY(&P2(rAbp,n,0), &P2(rAbp2,n,0));
	CONT(!(motSA.gTgl != 0));
	abpMotId[n] = abpMotId2[n];
  }

  free(chAbp2);
  free(rAbp2);
  free(del.l);
  free(abpIndL);
  if (motSA.gTgl != 0) { free(abpMotId2); }
}

// Elminate free actin segments. Because this code has neither polymerization
// nor depolymerization, they play no role here.
void DeleteFreeActin(void) {
  int n, k, side, cnt, *actInd, abpInd;

  Printf0("\n======================= Eliminating free actin segments "
			"========================\n");
  MALLOC(actInd, int, nAct);

  cnt = 0;
  FOR_ACT(n) {
	CONT(!(P2A(chAct,n,0,nChAc) > -1 || P2A(chAct,n,1,nChAc) > -1));
	actInd[n] = cnt;
	cnt++;
  }
  FOR_ACT(n) {
	for(k = 0; k < 2; k++) {
		CONT(!(P2A(chAct,n,k,nChAc) > -1));
		P2A(chAct,n,k,nChAc) = actInd[P2A(chAct,n,k,nChAc)]; 
	}
	for(k = 2; k < nChAc; k++) {
		CONT(!(P2A(chAct,n,k,nChAc) > -1));
		abpInd = P2A(chAct,n,k,nChAc);
		side = (n == P2A(chAbp,abpInd,0,nChAb)) ? 0 : 1;
		P2A(chAbp,abpInd,side,nChAb) = actInd[n];
	}
  }
  Printf0("%d free actin segments are eliminated.\n", nAct - cnt); 
  if (rheoWay > 0) {
	for(n = 0; n < appStraPar.c; n++) {
		appStraPar.l[n] = actInd[appStraPar.l[n]];
	}
	for(n = 0; n < meaStrePar.c; n++) {
		meaStrePar.l[n] = actInd[meaStrePar.l[n]];
	}
  }
  FOR_ACT(n) {
	CONT(!(P2A(chAct,n,0,nChAc) < 0 && P2A(chAct,n,1,nChAc) < 0));
	for(k = n; k < nAct - 1; k++) { 
		Copy1dArrayInt(&P2A(chAct,k,0,nChAc), &P2A(chAct,k + 1,0,nChAc), 
				nChAc);
		V3COPY(&P2(rAct,k,0), &P2(rAct,k + 1,0));
		actInd[k] = actInd[k + 1];
		fixAct[k] = fixAct[k + 1];
	}
	nAct--;
	n--;
  }
  free(actInd);
  MPI_Barrier(MPI_COMM_WORLD);
}

// Eliminate inactive ABPs if any. This is usally used to delete many inactive
// ABPs in a polymerized network. This should be executed before 
// ExtractConfig(), and it cannot be used at the mid of simulation.
void DeleteInactiveAbp(void) {
  int n, k, cnt, actInd, *chAbp2, side;
  double *rAbp2;

  Printf0("\n========================== Eliminating inactive ABPs "
			"===========================\n");

  MALLOC(chAbp2, int, nAbp*nChAb);
  MALLOC(rAbp2, double, nAbp*3);

  cnt = 0;
  // Information of active ABPs is stored in separate arrays.
  FOR_ABP(n) {
    if (P2A(chAbp,n,0,nChAb) > -1 && P2A(chAbp,n,1,nChAb) > -1) {
		Copy1dArrayInt(&P2A(chAbp2,cnt,0,nChAb), &P2A(chAbp,n,0,nChAb), nChAb);
		V3COPY(&P2(rAbp2,cnt,0), &P2(rAbp,n,0));
        for(k = 0; k < 2; k++) {
			actInd = P2A(chAbp,n,k,nChAb);
			side = FindAbpActinChain(actInd, n, 1);
			P2A(chAct,actInd,side,nChAc) = cnt;
        }
        cnt++;
    }
    else if (P2A(chAbp,n,0,nChAb) > -1 && P2A(chAbp,n,1,nChAb) < 0) {
		actInd = P2A(chAbp,n,0,nChAb);
		side = FindAbpActinChain(actInd, n, 1);
		P2A(chAct,actInd,side,nChAc) = -1;
    }
  }
  Printf0("%d of %d ABPs are deleted since they are in an inactive "
			"state!!\n", nAbp - cnt, nAbp);
  nAbp = cnt;
  nMot = 0;
  // Put it back to original array
  FOR_ABP(n) {
    Copy1dArrayInt(&P2A(chAbp,n,0,nChAb), &P2A(chAbp2,n,0,nChAb), nChAb);
    V3COPY(&P2(rAbp,n,0), &P2(rAbp2,n,0));
	if (P2A(chAbp,n,2,nChAb) == 2) { nMot++; }
  }
  free(chAbp2);
  free(rAbp2);
  MPI_Barrier(MPI_COMM_WORLD);
}

/*----------------------- Adding and deleting elements -----------------------*/

/*------------------- Severing filaments at the beginning --------------------*/

// Sever filaments whose length is greater than the minimum width of the 
// computational domain in order to preclude artifacts due to small domain. 
void SeverLongFilament(void) {
  int n, k, filaLen, curr, next, minDimDom = POS_LARGE_VALUE, ind;

  FOR_NDIM(k) { 
	CONT(dir2D == k);
    if ((int)dimDom[k] < minDimDom) { minDimDom = (int)dimDom[k]; }
  }
  FOR_ACT(n) {
	// Find a pointed end
    CONT(!(P2A(chAct,n,0,nChAc) > -1 && P2A(chAct,n,1,nChAc) < 0));
	filaLen = 2; 
	curr = P2A(chAct,n,0,nChAc);
	while (P2A(chAct,curr,0,nChAc) > -1) {
		curr = P2A(chAct,curr,0,nChAc);
		filaLen++;
		// If filament is longer than minimal width,the filament is severed.
		CONT(!(filaLen >= minDimDom));
		// Go backward by three segments
		for(k = 0; k < 3; k++) { curr = P2A(chAct,curr,1,nChAc); }
		next = P2A(chAct,curr,0,nChAc);
		// Sever it into two filaments
		P2A(chAct,curr,0,nChAc) = -1;
		P2A(chAct,next,1,nChAc) = -1;
		// Detach any ABP on the severed tip
		for(k = 0; k < 2; k++) {    
			ind = (k == 0) ? curr : next;
			DetachAbpOnActin(ind);
		}
		curr = next;
		filaLen = 1;
	}
  }
}

// Subroutine for SeverFilaOnPlane()
void SeverFilaOnPlaneSubroutine(int ind, int direc) {
  fixAct[ind]++;
  if (rheoWay > 0 && ((bulkRheoType == 0 && direc == dirNoPBC) 
		|| (bulkRheoType == 1 && direc == dirStr))) {
	fixAct[ind] += 10;
	InsertElement1dArrayWChk(meaStrePar.l,&meaStrePar.c,ind);
	InsertElement1dArrayWChk(appStraPar.l,&appStraPar.c,ind);
  }
}

// Sever all filaments crossing boundaries in a designated direction.
// After severing filaments, they are clamped to the boundaries.
// meaStrePar.l and appStraPar.l are filled here.
// mode = 1: clamped, 0: not clamped 
void SeverFilaOnPlane(int direc, int mode) {
  int n, k, l, ind, ind2, cnt, side, CS;
  double longLen;
  char namePlane[3];
 
  FOR_ACT(n) {
	ConvertRectDomainVector(&P2(rAct,n,0), 1);
  }
  FOR_ABP(n) {
	ConvertRectDomainVector(&P2(rAbp,n,0), 1);
  }
  switch(direc) {
    case 0: strcpy(namePlane, "yz"); break;
    case 1: strcpy(namePlane, "zx"); break;
    case 2: strcpy(namePlane, "xy"); break;
  }
  Printf0("\n==================== Severing filaments that cross %s plane "
			"====================\n", namePlane);
  cnt = 0;

  // Delete all ABP bonds near boundaries.
  longLen = 1.0;
  for(n = 0; n < NK_ABP; n++) {
	if (abpF.len[n].n + 0.5 * actF.dia > longLen) {
		longLen = abpF.len[n].n + 0.5 * actF.dia;
	}
  }
  FOR_ABP(n) {
    CONT(!(P2(rAbp,n,direc) >= rGrid[direc][nGrid[direc] - 1] - longLen
	        || P2(rAbp,n,direc) < rGrid[direc][0] + longLen));
	cnt++;
	DetachActinOnAbp(n);
  }
  Printf0("%d ABPs located near %s plane are unbound.\n", cnt, namePlane);

  // Fill in meaStrePar.l and appStraPar.l
  cnt = 0;
  if (rheoWay > 0 && ((bulkRheoType == 0 && direc == dirNoPBC)
        || (bulkRheoType == 1 && direc == dirStr))) {
	meaStrePar.c = 0;
	appStraPar.c = 0;
  }
  // Find filaments crosslinking boundaries and sever/clamp them.
  FOR_ACT(n) {
	CONT(!(P2(rAct,n,direc) >= rGrid[direc][nGrid[direc] - 1] - 1.0 
			|| P2(rAct,n,direc) <= rGrid[direc][0] + 1.0 ));
	CS = (P2(rAct,n,direc) > rGrid[direc][0] + dimDomH[direc]) ? 1 : 0;
	for(k = 0; k < 2; k++) {
		ind = P2A(chAct,n,k,nChAc);
		CONT(!(ind > -1));
		CONT(!(fabs(P2(rAct,n,direc) - P2(rAct,ind,direc)) > dimDomH[direc]));
		if (P2A(chAct,ind,k,nChAc) > -1
				&& P2A(chAct,n,1 - k,nChAc) > -1) {
			if (mode == 1) {
				for(l = 0; l < 2; l++) {
					ind2 = (l == 0) ? n : ind;
					side = (l == 0) ? 1 - k : k;		
					CONT(!(fixAct[P2A(chAct,ind2,side,nChAc)] < 0));
					fixAct[ind2] = direc * 2 + 1;
					if ((l == 0 && CS == 1) || (l == 1 && CS == 0)) {
						SeverFilaOnPlaneSubroutine(ind2, direc);
					}
				}
			}
		}
		else {
			if (rheoWay > 0) { 
				for(l = 0; l < 2; l++) {
					ind2 = (l == 0) ? n : ind;
					DeleteElement1dArray(meaStrePar.l, &meaStrePar.c, ind2);
					DeleteElement1dArray(appStraPar.l, &appStraPar.c, ind2);
				}
			}
		}
		cnt++;
		DetachAbpOnActin(n);
		DetachAbpOnActin(ind);
		P2A(chAct,n,k,nChAc) = -1;
		P2A(chAct,ind,1 - k,nChAc) = -1;
	}
  }
  Printf0("%d filaments crossing %s plane are severed.\n",
			 cnt, namePlane);
  if (rheoWay > 0 && ((bulkRheoType == 0 && direc == dirNoPBC)
        || (bulkRheoType == 1 && direc == dirStr))) {
	qsort(appStraPar.l, appStraPar.c, sizeof(int), CompInt);
	qsort(meaStrePar.l, meaStrePar.c, sizeof(int), CompInt);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

/*------------------- Severing filaments at the beginning --------------------*/

/*----------------- Breaking chains between actins and ABPs ------------------*/

// Detach all ABP bonds from the actin segment
void DetachAbpOnActin(int actInd) {
  int n, abpInd;

  for(n = 2; n < nChAc; n++) {
	CONT(!(P2A(chAct,actInd,n,nChAc) > -1));
	abpInd = P2A(chAct,actInd,n,nChAc);
	if (P2A(chAbp,abpInd,0,nChAb) == actInd) {
		P2A(chAbp,abpInd,0,nChAb) = P2A(chAbp,abpInd,1,nChAb);
	}
	P2A(chAbp,abpInd,1,nChAb) = -1;
	P2A(chAct,actInd,n,nChAc) = -1;
  }
}

// Detach all actin bonds from the ABP
void DetachActinOnAbp(int abpInd) {
  int n, actInd, side;

  for(n = 0; n < 2; n++) {
	CONT(!(P2A(chAbp,abpInd,n,nChAb) > -1));
	actInd = P2A(chAbp,abpInd,n,nChAb);
	side = FindAbpActinChain(actInd, abpInd, 1);
	P2A(chAct,actInd,side,nChAc) = -1;
	P2A(chAbp,abpInd,n,nChAb) = -1;
  }
}

/*----------------- Breaking chains between actins and ABPs ------------------*/

/*------------------------ Inspect chain information -------------------------*/

// Check the integrity of loaded network configuration in terms of 
// chain information.
int InspectChainList(void) {
  int n, k, l, CS, CS2, ind1, ind2, ind3;

  CS = 1;
  FOR_ACTME(n) {
	// Check chain information between actin segments in "act.ch"
	ind1 = act.id[n];
	for(k = 0; k < 2; k++) {
		CONT(!(P2A(act.ch,n,k,nChAc) > -1));
		ind2 = P2A(act.ch,n,k,nChAc);
		ind3 = iAct[ind2];
		CONT(!(P2A(act.ch,ind3,1 - k,nChAc) != ind1));
		Printf("Error: chain information between actin segment %d and"
				" %d is wrong at rank %d.\n", ind1, ind2, rank);
		Printf("%d: ", ind1);
		for(l = 0; l < nChAc; l++) {
			Printf("%d\t", P2A(act.ch,n,l,nChAc));
		}
		Printf("\n%d: ", ind2);
		for(l = 0; l < nChAc; l++) {
			Printf("%d\t", P2A(act.ch,ind3,l,nChAc));
		}
		Printf("\n\n");
		CS = 0;
	}
	// Check chain information between actin segment and ABP in "act.ch"
	for(k = 2; k < nChAc; k++) {
		CONT(!(P2A(act.ch,n,k,nChAc) > -1));
		ind2 = P2A(act.ch,n,k,nChAc);
		ind3 = iAbp[ind2];
		CONT(!(P2A(abp.ch,ind3,0,nChAb) != ind1 
				&& P2A(abp.ch,ind3,1,nChAb) != ind1));
		Printf("Error: chain information between actin segment %d and"
				" ABP %d is wrong at rank %d.\n", ind1, ind2, rank);
		Printf("%d: ", ind1);
		for(l = 0; l < nChAc; l++) {
			Printf("%d\t", P2A(act.ch,n,l,nChAc));
		}
		Printf("\n%d: ", ind2);
		for(l = 0; l < nChAb; l++) {
			Printf("%d\t", P2A(abp.ch,ind3,l,nChAb));
		}
		Printf("\n\n");
		P2A(act.ch,n,k,nChAc) = -1;
		CS = 0;
	}
  }
  // Check chain information between actin segment and ABP in "abp.ch"
  FOR_ABPME(n) {
	ind1 = abp.id[n];
	for(k = 0; k < 2; k++) {
		CONT(!(P2A(abp.ch,n,k,nChAb) > -1));
		ind2 = P2A(abp.ch,n,k,nChAb);
		ind3 = iAct[ind2];
		CS2 = FindAbpActinChain(ind3, ind1, 0);
		CONT(!(CS2 == -1));
		Printf("Error: chain information between actin segment %d and"
				" ABP %d is wrong at rank %d.\n", ind2, ind1, rank);
		Printf("%d: ", ind2);
		for(l = 0; l < nChAc; l++) {
			Printf("%d\t", P2A(act.ch,ind3,l,nChAc));
		}
		Printf("\n%d: ", ind1);
		for(l = 0; l < nChAb; l++) {
			Printf("%d\t", P2A(abp.ch,n,l,nChAb));
		}
		Printf("\n\n");
		CS = 0;
	}
  }
  return CS;
}

// Check the integrity of loaded network configuration in terms of chain length
int InspectChainLength(void) {
  int n, k, ind1, ind2, ind3, CS = 1;
  double len, dist;

  // Check chain lengths using "act.ch"
  FOR_ACTME(n) {
	ind1 = act.id[n];	
	for(k = 0; k < nChAc; k++) {
		ind2 = P2A(act.ch,n,k,nChAc);
		CONT(ind2 < 0);
		if (k < 2) { 
			ind3 = iAct[ind2];
			len = CalcDist(&P2(act.r,n,0), &P2(act.r,ind3,0), 0); 
		}
		else { 
			ind3 = iAbp[ind2];
			len = CalcDistActinAbp(ind1, ind2, 0);
		}
		if (k < 2) { dist = AVG2(1., 1.); }
		else { dist = 0.5 * actF.dia + abpF.len[K_ABP(ind3)].n; }

   		CONT(!(len < dist - 0.5 || len > dist + 0.5));
		if (k < 2) {
			Printf("Error: distance between actin segment %d and %d is "
					"wrong at rank %d: %lf (%lf)\n", 
					ind1, ind2, rank, len, dist);
			Printf("%d: %lf\t%lf\t%lf\n", ind1, P2(act.r,n,0), 
					P2(act.r,n,1), P2(act.r,n,2));
			Printf("%d: %lf\t%lf\t%lf\n", ind2, P2(act.r,ind3,0), 
					P2(act.r,ind3,1), P2(act.r,ind3,2));
		}
		else {
			Printf("Error: distance between actin segment %d and ABP"
					" %d is wrong at rank %d: %lf (%lf)\n", 
					ind1, ind2, rank, len, dist);
			Printf("%d: %lf\t%lf\t%lf\n", ind1, P2(act.r,n,0), 
					P2(act.r,n,1),	P2(act.r,n,2));
			Printf("%d: %lf\t%lf\t%lf\n", ind2, P2(abp.r,ind3,0), 
					P2(abp.r,ind3,1), P2(abp.r,ind3,2));
		}
		CS = 0;
	}
  }
  // Check chain lengths using "abp.ch"
  FOR_ABPME(n) {
	ind1 = abp.id[n];
	for(k = 0; k < 2; k ++) {
		ind2 = P2A(abp.ch,n,k,nChAb);
		CONT(ind2 < 0);
		ind3 = iAct[ind2];
		len = CalcDistActinAbp(ind2, ind1, 0);
		dist = 0.5 * actF.dia + abpF.len[K_ABP(n)].n;
		CONT(!(len < dist - 0.5 || len > dist + 0.5));
		Printf("Error: distance between actin segment %d and ABP"
				" %d is wrong at %d: %lf (%lf)\n", ind2, ind1, rank, len, dist);
		Printf("%d: %lf\t%lf\t%lf\n", ind2, P2(act.r,ind3,0), 
				P2(act.r,ind3,1), P2(act.r,ind3,2));
		Printf("%d: %lf\t%lf\t%lf\n", ind1, P2(abp.r,n,0), 
				P2(abp.r,n,1), P2(abp.r,n,2));
		CS = 0;
	}
  }
  return CS;
}

// Inspect chain length and information by calling InspectChainList() and
// InspectChainLength().
void InspectChainListAndLength(void) {
  int CS1, CS2;
  Printf0("\n================ Inspecting the integrity of chain information "
		"=================\n");
  CS1 = InspectChainList();
  CS2 = InspectChainLength();
  if (CS1 != 1) { Printf("Chain list has a problem!!\n"); }
  if (CS2 != 1) { Printf("Chain length has a problem!!\n"); }
  MPI_Barrier(MPI_COMM_WORLD);
}

int BinaryPackActinChainArray(int actInd) {
  int k, ch, cntAbp;

  ch = 0;
  cntAbp = 0;
  for(k = 0; k < nChAc; k++) {
	CONT(!(P2A(act.ch,iAct[actInd],k,nChAc) > -1));
	if (k < 2) { ch += (int)(pow(2, k + 1)); }
	else { cntAbp++; }
  }
  k = 3;
  while(cntAbp != 0) {
	if (cntAbp % 2 == 1) { ch += (int)(pow(2, k)); }
	cntAbp /= 2;
	k++;
  }
  if (gTglActCapAll != 0) { 
	if (act.cap[iAct[actInd]] > -1) {
		ch += 1;
	}
  }
  return ch;
}

/*------------------------ Inspect chain information -------------------------*/

/*--------------------------- Pack information -------------------------------*/

// 0: actin0, 2: stalling at barbed end, 4: stalling by traffic jam; 
// 6: stalling by forces
// 1, 3, 5, 7: same for actin 1
// 8: the kind of ABP
int BinaryPackAbpChainArray(int abpInd) {
  int k, ch, side, actInd, locActInd, locAbpInd, locNextActInd;
  double pMotW, critPMotW;

  critPMotW = 1.0e-10;
  ch = 0;
  locAbpInd = iAbp[abpInd];
  for(k = 0; k < 2; k++) {
	actInd = P2A(abp.ch,locAbpInd,k,nChAb);
	// Check whether any actin is bound on ABP
	if (actInd > -1) { 
		ch += (int)(pow(2, k)); 
		locActInd = iAct[actInd];
		CONT(locActInd < 0);
	}
	else { continue; }
	locNextActInd = iAct[P2A(act.ch,locActInd,0,nChAc)];
	// Check traffic jam
	side = FindAbpActinChain(locActInd, abpInd, 0);
	side = (side - 2) / nChAcY;
	if (side < nChAcX - 1) {
		side = (((side - 2) / nChAcY) + 1) * nChAcY + 2;
		side = FindElementArray(&P2A(act.ch,locActInd,side,nChAc), nChAcY, 
				-1, 0, 1);
	}
	else {
		if (locNextActInd > -1) {
			side = FindElementArray(&P2A(act.ch,locNextActInd,2,nChAc), nChAcY, 
					-1, 0, 1);
			// Check whether ABP is located at the barbed ends
			if (P2A(act.ch,locNextActInd,0,nChAc) < 0) {
				ch += (int)pow(2, k + 2);
			}
		}
	}
	if (side < 0) { 
		ch += (int)pow(2., k + 4);
	}
	CONT(K_ABP(locAbpInd) != 2 || motWalk.gTgl == 0);
	// Check stalling by forces
	pMotW = UpdateMotorWalkingSubroutine(abpInd, k);
	if (pMotW < critPMotW) {
		ch += (int)pow(2., k + 6);
	}
  }
  ch += (int)pow(2., K_ABP(locAbpInd) + 8);
  if (ISMTF(K_ABP(locAbpInd))) {
	if (abp.mId[locAbpInd] > -1) {
		ch += (int)pow(2., abp.mId[locAbpInd] + 11);
	}
  }
  return ch;
}

/*--------------------------- Pack information -------------------------------*/

/*---------------- Find supportive and percolated structures -----------------*/

void FindSupportiveFramework(void) {
  int mm, m, n, k, curr, chkPnt;
  int *nActAll, *nAbpAll, *hiFabp;
  double *arr, *recFabpAll;

  MALLOC(nActAll,int,nCpu);
  MALLOC(nAbpAll,int,nCpu);
  if (rank == 0) {
    MALLOC(chAct,int,nAct * nChAc);
    MALLOC(chAbp,int,nAbp * nChAb);
	MALLOC(recFabpAll, double, nAbp);
	MALLOC(fixAct,int,nAct);
  }
  MPI_Gather(&nActMe, 1, MPI_INT, nActAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  Gather2dArrayInt(nActMe, nActAll, nChAc, act.id, act.ch, chAct);
  Gather2dArrayInt(nActMe, nActAll, 1, act.id, act.fix, fixAct);
  MPI_Gather(&nAbpMe, 1, MPI_INT, nAbpAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  Gather2dArrayInt(nAbpMe, nAbpAll, nChAb, abp.id, abp.ch, chAbp);
  Gather2dArrayDouble(nAbpMe, nAbpAll, 1, abp.id, 
		((kindFindSupp == 0) ? recAbp.bendF : recAbp.sprF), recFabpAll);

  if (rank == 0) {
	MALLOC(hiFabp,int,nAbp);
	MALLOC(arr,double,nAbp*2);
	memset(hiFabp, -1, sizeof(int) * nAbp);
	memset(iSupp, -1, sizeof(int) * (nAct + nAbp));
	for(n = 0; n < nAbp; n++) { 
		if (P2A(chAbp,n,0,nChAb) < 0 || P2A(chAbp,n,1,nChAb) < 0) {
			recFabpAll[n] = 0.;
		}
		V2SET(&P2A(arr,n,0,2), recFabpAll[n], (double)n);
	}
	qsort(arr, nAbp, 2 * sizeof(double), CompDbl);
	for(n = 0; n < (int)(porFindSupp * (double)nAbp); n++) {
		hiFabp[(int)P2A(arr,nAbp - n - 1,1,2)] = 1;
	}
	for(m = 0; m < nAbp; m++) {
		CONT(hiFabp[m] < 0);
		iSupp[nAct + m] = 1;
		for(n = 0; n < 2; n++) {
			CONT(P2A(chAbp,m,n,nChAb) < 0);
			for(k = 0; k < 2; k++) {
				curr = P2A(chAbp,m,n,nChAb);
				while(curr > -1) {
					for(mm = 0; mm < nChAc - 2; mm++) {
						CONT(P2A(chAct,curr,mm + 2,nChAc) < 0);
						if (hiFabp[P2A(chAct,curr,mm + 2,nChAc)] > -1) {
							chkPnt = curr;
							break;
						}
					}
					if (fixAct[curr] > -1) {
						chkPnt = curr;
					}
					curr = P2A(chAct,curr,k,nChAc);
				}
				if (chkPnt != P2A(chAbp,m,n,nChAb)) {
					curr = chkPnt;
					do {
						iSupp[curr] = 1;
						curr = P2A(chAct,curr,1 - k,nChAc);
					} while(curr != P2A(chAbp,m,n,nChAb));
				}
				iSupp[chkPnt] = 1;
			}
		}
	}	
	free(arr);
	free(hiFabp);
	free(chAct);
	free(chAbp);
	free(fixAct);
	free(recFabpAll);
  }
  free(nActAll);
  free(nAbpAll);
  MPI_Bcast(iSupp, nAct + nAbp, MPI_INT, 0, MPI_COMM_WORLD);
}

void FindPercolatedFilamentsSubroutine(ListInt *filaPath, ListInt *percPath, 
		ListInt *dnF, ListInt2 *recPath, int *percPathLen, int *iFila, 
		int *chkIfila, int nFila, int direc, int *maxCol) {

  int n, k, l, CS;
  int curr, start, *chkL, thres, abpInd, side, actInd;
  ListInt segPath;

  // If the pathway is longer than [1.8 x domain dimension], search is ended.
  thres = (int)(dimDom[(bulkRheoType == 0) ? dirNoPBC : direc] * 1.2);
  MALLOC(chkL, int, nFila);
  MALLOC(segPath.l, int, nAct);
  memset(chkL, -1, sizeof(int) * nFila);
  start = P2A(filaPath->l,filaPath->c - 1,1,2);
  // 0: barbed-end direction, 1: pointed-end direction
  for (n = 0; n < 2; n++) {
	curr = start;
	segPath.c = 1;
	segPath.l[0] = curr;
	// segPath: remember all actin indices for one segment
	// percPath: store all actin indices for pathways
    if (n == 1) { 
		curr = P2A(chAct,curr,n,nChAc); 
        segPath.l[segPath.c] = curr;
        (segPath.c)++;
	}
    while(curr > -1 && segPath.c + percPath->c < thres) {
        for(k = 0; k < nChAc - 2; k++) {
            abpInd = P2A(chAct,curr,k + 2,nChAc);
            CONT(abpInd < 0);
			// Connections made by self-assembled motors are ignored.
			CONT(ISMTF(P2A(chAbp,abpInd,2,nChAb)));
			// Ignore ABP that already passed.
			CONT(nAct + abpInd == percPath->l[percPath->c - 1]);
			// side: which ABP arm is connected
            side = (P2A(chAbp,abpInd,0,nChAb) == curr) ? 0 : 1;
			// actInd: the actin index on the other arm
            actInd = P2A(chAbp,abpInd,1 - side, nChAb); 
            CONT(actInd < 0);
			// Ignore a filament that already passed
            CONT(chkL[iFila[actInd]] > -1);
			// If a new filament was already counted, it is ignored.
            for(l = 0; l < filaPath->c; l++) {
                BREAK(P2A(filaPath->l,l,0,2) == iFila[actInd]);
            }
            CONT(l < filaPath->c);
			// chkL: tag for the filaments
			chkL[iFila[actInd]] = 1;
			// filaPath: remember actin filament indices that passed
            V2SET(&P2A(filaPath->l,filaPath->c,0,2), iFila[actInd], actInd);
			(filaPath->c)++;

            for(l = 0; l < segPath.c; l++) {
                percPath->l[percPath->c + l] = segPath.l[l];
            }   
            (percPath->c) += segPath.c;
		    percPath->l[percPath->c] = abpInd + nAct;
            (percPath->c) += 1;
            FindPercolatedFilamentsSubroutine(filaPath, percPath, dnF, recPath,
                    percPathLen, iFila, chkIfila, nFila, direc, maxCol); 
			// Follow the previous pathways backwards
			while(percPath->l[percPath->c] != start) {
				percPath->l[percPath->c] = -1;
				(percPath->c)--;
			}
			percPath->l[percPath->c] = -1;
			(filaPath->c)--;
        }
		if (P2A(chAct,curr,n,nChAc) < 0) {
			// Compare whether they are ending points
			CS = FindElementArray(dnF->l, dnF->c, curr, 1, 2);

			if (CS > -1) {
				for(k = 0; k < percPath->c; k++) {
					iPerc[percPath->l[k]] = 1;
				}
				for(k = 0; k < segPath.c; k++) {
					iPerc[segPath.l[k]] = 1;
				}
				// recPath: store filament indices
				MALLOC(recPath->l[recPath->c], int, nFila + 1);
				memset(recPath->l[recPath->c], -1, sizeof(int) * (nFila + 1));
				recPath->l[recPath->c][0] = percPath->c + segPath.c;
				percPathLen[percPath->c + segPath.c]++;
				// store indices of filaPath in recPath.
				for(k = 0; k < filaPath->c; k++) {
					recPath->l[recPath->c][k + 1] 
							= chkIfila[P2A(filaPath->l,k,0,2)];
				}
				if (filaPath->c > *maxCol) { *maxCol = filaPath->c; }
				(recPath->c)++;
			}
		} 
		// Jump to a next actin point
        curr = P2A(chAct,curr,n,nChAc);
		// segPath: store actin point indices
        segPath.l[segPath.c] = curr;
        (segPath.c)++;
    }
  }
  free(chkL);
  free(segPath.l);
}

void FindPercolatedFilaments(int direc, int period) {
  int n, k, upFLind, CS, nFila, sizUpFL;
  int *nActAll, *nAbpAll, *iFila, *chkIfila;
  int maxCol, cntArr[2], *cntArrAll, *iPercAll, sumRecPathC, tag = 0;
  int *percPathLenAll, *percPathLen, maxPathLen, minPathLen;
  ListInt upF, dnF, percPath, filaPath;
  ListInt2 recPath;
  FILE *fOut;

  if (pbc[direc] != 0) {
	Printf0("Error: Periodic boundary condition exists in the direction!\n\n");
	exit(-1);
  }

  MALLOC(nActAll,int,nCpu);
  MALLOC(nAbpAll,int,nCpu);
  MALLOC(chAct,int,nAct * nChAc);
  MALLOC(chAbp,int,nAbp * nChAb);
  MALLOC(iFila,int,nAct);
  MALLOC(chkIfila,int,nAct);
  MALLOC(fixAct,int,nAct);
  MALLOC(percPath.l,int,nAct); 
  MALLOC(filaPath.l,int,nAct * 2); 
  MALLOC(percPathLen,int,nAct);
  if (rank == 0) { 
	MALLOC(percPathLenAll, int, nAct * nCpu);
	MALLOC(iPercAll, int, (nAct + nAbp) * nCpu);
	MALLOC(cntArrAll, int, nCpu * 2);
  }
  maxPathLen = NEG_LARGE_VALUE;
  minPathLen = POS_LARGE_VALUE;
  maxCol = NEG_LARGE_VALUE;
  memset(percPathLen, 0, sizeof(int) * nAct);
  memset(iPerc, 0, sizeof(int) * (nAct + nAbp));
  recPath.c = 0;

  MPI_Gather(&nActMe, 1, MPI_INT, nActAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  Gather2dArrayInt(nActMe, nActAll, nChAc, act.id, act.ch, chAct);
  Gather2dArrayInt(nActMe, nActAll, 1, act.id, act.iF, iFila);
  Gather2dArrayInt(nActMe, nActAll, 1, act.id, act.fix, fixAct);
  MPI_Gather(&nAbpMe, 1, MPI_INT, nAbpAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  Gather2dArrayInt(nAbpMe, nAbpAll, nChAb, abp.id, abp.ch, chAbp);

  // chkIfila has the original index of filaments
  // iFila has the resorted index
  RecordConnectivitySubroutine1(chkIfila, iFila, &nFila);

  // Broadcast chAct, chAbp, fixAct, and iFila to all
  MPI_Bcast(chAct, nAct * nChAc, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(chAbp, nAbp * nChAb, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(fixAct, nAct, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(iFila, nAct, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(chkIfila, nAct, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&dimDom[direc], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // upF.l: list of filaments attached to the top boundary
  // dnF.l: list of filaments attached to the bottom boundary
  MALLOC(upF.l, int, nFila * 2);  
  MALLOC(dnF.l, int, nFila * 2);  
  upF.c = 0;
  dnF.c = 0;
  // Fill in "upF.l" and "dnF.l"
  FOR_ACT(n) {
	CONT((P2A(chAct,n,0,nChAc) > -1 && P2A(chAct,n,1,nChAc) > -1)
			|| (P2A(chAct,n,0,nChAc) < 0 && P2A(chAct,n,1,nChAc) < 0));
	if (fixAct[n] % 10 == 2 * direc + 2) {
		CS = FindElementArray(upF.l, upF.c, iFila[n], 0, 2);
		if (CS == -1) { 
			V2SET(&P2A(upF.l,upF.c,0,2), iFila[n], n);
			(upF.c)++;
		}
	}
	if (fixAct[n] % 10 == 2 * direc + 1) { 
		CS = FindElementArray(dnF.l, dnF.c, iFila[n], 0, 2);
		if (CS == -1) { 
			V2SET(&P2A(dnF.l,dnF.c,0,2), iFila[n], n);
			(dnF.c)++;
		}
	}
  }
  qsort(upF.l, upF.c, 2 * sizeof(int), CompInt);
  qsort(dnF.l, dnF.c, 2 * sizeof(int), CompInt);
//////////////////////////////////////////////////////////////
  MALLOC2(recPath.l,int,nFila * upF.c * 100);
//////////////////////////////////////////////////////////////
  sizUpFL = (upF.c / nCpu) + (rank < (upF.c % nCpu) ? 1 : 0);
  upFLind = 0; 
  for(n = 0; n < rank; n++) { 
  	upFLind += (upF.c / nCpu) + (n < (upF.c % nCpu) ? 1 : 0);
  }
  for(n = 0; n < sizUpFL; n++) {
	filaPath.c = 1;
	V2COPY(&P2A(filaPath.l,0,0,2), &P2A(upF.l,n + upFLind,0,2));
	percPath.c = 0;
	FindPercolatedFilamentsSubroutine(&filaPath, &percPath, &dnF, &recPath, 
			percPathLen, iFila, chkIfila, nFila, direc, &maxCol);
	Check2dArraySize(&recPath, &recPath.siz, nFila + 1, 1);
  }

  V2SET(cntArr, recPath.c, maxCol);
  MPI_Gather(percPathLen, nAct, MPI_INT, percPathLenAll, nAct, MPI_INT, 0, 
		MPI_COMM_WORLD);
  MPI_Gather(iPerc, nAct + nAbp, MPI_INT, iPercAll, nAct + nAbp, 
		MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(cntArr, 2, MPI_INT, cntArrAll, 2, MPI_INT, 0, MPI_COMM_WORLD);
  if (rank == 0) {
	// iPerc: whether an element belongs to percolating pathways or not
	for(n = 0; n < nAct + nAbp; n++) {
		for(k = 1; k < nCpu; k++) {
			if (iPercAll[k * (nAct + nAbp) + n] > iPerc[n]) {
				iPerc[n] = iPercAll[k * (nAct + nAbp) + n];
			}
		}
	}
	// percPathLen: the number of pathways for each length
	for(n = 0; n < nAct; n++) {
		for(k = 1; k < nCpu; k++) {
			percPathLen[n] += percPathLenAll[k * nAct + n];
		}
	}
	// sumRecPathC: the total number of athways
	sumRecPathC = recPath.c;
	for(n = 1; n < nCpu; n++) {
		if (P2A(cntArrAll,n,1,2) > maxCol) { 
			maxCol = P2A(cntArrAll,n,1,2);
		}	
		sumRecPathC += P2A(cntArrAll,n,0,2);
	}
	for(n = 0; n < nAct; n++) {
		CONT(percPathLen[n] == 0);
		if (maxPathLen < n) { maxPathLen = n; }
		if (minPathLen > n) { minPathLen = n; }
	}
	fOut = fopen(GenFileName("PercLen"), "a");
	// Record the number of lines and the instantaneous width of domain.
	fprintf(fOut, "%d\t%g\t0\n", ((maxPathLen + 1 >  minPathLen) ? 
			maxPathLen - minPathLen + 1 : 0),
			dimDom[(bulkRheoType == 0) ? dirNoPBC : dirStr]);
	// Record percolation information
	for(n = minPathLen; n < maxPathLen + 1; n++) {
		fprintf(fOut, "%d\t%d\t%d\n", n - minPathLen, n, percPathLen[n]);
	}
	fclose(fOut);
  }
  MPI_Bcast(iPerc, nAct + nAbp, MPI_INT, 0, MPI_COMM_WORLD);

  // By using MPI_Send and MPI_Recv, each CPU records data in PercPath 
  // sequentially, and they can also synchronize maxCol.
  if (rank != 0) {	
	MPI_Recv(cntArr, 2, MPI_INT, rank - 1, tag, MPI_COMM_WORLD, &status);
	maxCol = cntArr[1];
  }
  else { cntArr[0] = 0; }

  fOut = fopen(GenFileName("PercFila"), "a");
  if (rank == 0) {
	fprintf(fOut, "%d\t%g\t%d\t%d\t", sumRecPathC, dimDom[direc], 
			nFila, ((maxCol < 0) ? 0 : maxCol));
	Fprintf1dFillerInt(fOut, 0, maxCol - 2, 0);
  }
  for(n = 0; n < recPath.c; n++) {
	fprintf(fOut, "%d\t", n + cntArr[0]);
	Fprintf1dArrayInt(fOut, recPath.l[n], maxCol + 1, 0);
  }
  fclose(fOut);
  if (rank < nCpu - 1) { 
	V2SET(cntArr, recPath.c + cntArr[0], maxCol);
	MPI_Send(cntArr, 2, MPI_INT, rank + 1, tag, MPI_COMM_WORLD);
  }

  for(n = 0; n < recPath.c; n++) {
	free(recPath.l[n]);
  }

  if (rank == 0) { 
	free(percPathLenAll);
	free(iPercAll);
	free(cntArrAll);
  }
  free(upF.l);  
  free(dnF.l);  
  free(recPath.l);
  free(nActAll);
  free(nAbpAll);
  free(chAct);
  free(chAbp);
  free(iFila);
  free(chkIfila);
  free(fixAct);
  free(percPath.l);
  free(filaPath.l); 
  free(percPathLen);
}
