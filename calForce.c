// ##################################################
// #   calForce.c - finally revised on Dec 2018     #
// #   coded by Taeyoon Kim                         #
// #   Copyright (C) 2005 - 2018, Taeyoon Kim,      #
// #   All rights reserved.                         #
// ##################################################
// This file contains functions which calculate forces.

/*----------------------- Related to thermal fluctuation ---------------------*/

inline void CalcBrownForcesSubroutine(double *fBr, double maxFbr, 
		int ind1, int ind2, double *rand) {

  if (ind1 % 2 == 0) {
   	genrand_gauss(&rand[0], &rand[1]);
    P2(fBr,ind2,0) = rand[0] * maxFbr;
    P2(fBr,ind2,1) = rand[1] * maxFbr;
   	genrand_gauss(&rand[0], &rand[1]);
    P2(fBr,ind2,2) = rand[0] * maxFbr;
  }
  else {
    P2(fBr,ind2,0) = rand[1] * maxFbr;
   	genrand_gauss(&rand[0], &rand[1]);
    P2(fBr,ind2,1) = rand[0] * maxFbr;
    P2(fBr,ind2,2) = rand[1] * maxFbr;
  }
}
// Apply Brownian forces (thermal fluctuation) to actin, ACP, and motor.
void CalcBrownForces(void) {
  int n, k, *pArr, kind, cnt;
  double rand[2], maxFbr;

  if (gTglActTherm != 0) {
	cnt = 0;
	FOR_ACTME(n) {
		CONT(ISACTM(n));
		CalcBrownForcesSubroutine(act.fBr, act.maxFbr[0], cnt, n, rand);
		cnt++;
	}
  }
  if ((gTglAcpTherm != 0 && nAbp - nMot > 0) 
		|| (gTglMotTherm != 0 && nMot > 0)) {
	cnt = 0;
	FOR_ABPME(n) {
		CONT(ISABPIM(n));
		kind = K_ABP(n);
		if (kind == 2) { CONT(gTglMotTherm == 0); }
		else { CONT(gTglAcpTherm == 0); }
		CalcBrownForcesSubroutine(abp.fBr, abp.maxFbr[kind], cnt, n, rand);
		cnt++;
	}
  }
}

/*----------------------- Related to thermal fluctuation ---------------------*/

/*------------------------ Related to repulsive forces -----------------------*/

void CalcRepulsiveForcesSubSubSubroutine(double rPnt[][NDIM], 
		int begin, int end) {
  int n, k;

  for(n = begin; n <= end; n++) {
    FOR_NDIM(k) {
        CONT(pbc[k] != 1);
        CONT(!(fabs(rPnt[n][k] - rPnt[0][k]) > dimDomH[k]));
        rPnt[n][k] += dimDom[k]
                * ((rPnt[n][k] < rPnt[0][k]) ? 1. : -1.);
	}
  }
}

double CalcRepulsiveForcesSubSubroutine(double rPnt[][NDIM], double *dr, 
		double *ratio, double dist, int mode) {
  int CS, dim;
  double len;

  dim = 2; 
  if (mode == 0) {
	CalcRepulsiveForcesSubSubSubroutine(rPnt, 1, 3);
	len = CalcSegSegDist(&rPnt[0], &rPnt[2], dr, ratio, 0);
  }
  else if (mode == 1) {
	CalcRepulsiveForcesSubSubSubroutine(rPnt, 1, 2);
	len = CalcSegPntDist(&rPnt[0], rPnt[2], dr, ratio);
  }
  else if (mode == 2) {
	CalcRepulsiveForcesSubSubSubroutine(rPnt, 2, 2);
	len = CalcVecDist(dr, rPnt[0], rPnt[2], 0);
  }
  return len;
}
// mode = 0: actin and actin
//        1: actin and ABP
//        2: ABP and ABP
void CalcRepulsiveForcesSubroutine(int *ind, int *nlPnt, double rPnt[][NDIM], 
		double *f_p[], double *dist, int mode) {
  int n, k, chkErr, dim;
  double f, len, dist2, dot, fi[NDIM], dr[NDIM], normDir[NDIM], ratio[4];
  double stf, rel, *pR;

  dim = 2; 
  dist2 = AVG2(dist[0], dist[1]);
  len = CalcRepulsiveForcesSubSubroutine(rPnt, dr, ratio, dist2, mode);
  if (mode == 4 && actDgd.gTgl != 0) {
	UpdateActinDegradation(ind, len, dim);
  }
  if (len == 0.) { return; }

  // If the distance between two segments is small than the sum of
  // their radii (overlapping), repulsive forces should act.
  if (len < dist2) {
	if (mode == 0) { stf = actF.rep.stf; }
	else if (mode == 1 || mode == 2) { stf = abpF.rep.stf; }
	else { stf = stfRepMbNuc; }

	if (mode < 3) {
		f = stf * (dist2 - len);
	}
	chkErr = CheckLargeForce(f, 1);

	if (chkErr != 1) {
		RecordErrorRepForce(ind, nlPnt, rPnt, f, len, ratio, mode);
	}
	f /= len;
	VS3COPY(fi, dr, f);
	if (mode == 0) {
		VVS3ADD(f_p[0], fi, REVSIGN(1. - ratio[0]));
		VVS3ADD(f_p[1], fi, REVSIGN(ratio[0]));
		VVS3ADD(f_p[2], fi, 1. - ratio[1]);
		VVS3ADD(f_p[3], fi, ratio[1]);
	}
	else if (mode == 1) {
		VVS3ADD(f_p[0], fi, REVSIGN(1. - ratio[0]));
		VVS3ADD(f_p[1], fi, REVSIGN(ratio[0]));
		VV3ADD(f_p[2], fi);
	}
	else if (mode == 2) {
		VV3SUB(f_p[0], fi);
		VV3ADD(f_p[2], fi);
	}
  }
}

// Calculate volume-exclusion effects between actin segments.
void CalcRepulsiveForces(void) {
  int n, k, *nlPnt, ind[4], ind2[2], mode, abpInd, abpInd2;
  double rPnt[4][NDIM], *r_p[4], *f_p[4], dist[2];
double len, dist2, f, fi[NDIM], ratio[2], dr[NDIM];

  nlPnt = neigh.l;
  for(n = 0; n < neigh.c; n++) {
	if (n > 0) { nlPnt += 2; }
	// From the neighboring list, required information is prepared
	CONT(actF.rep.facStf == 0. && nlPnt[0] < nAct && nlPnt[1] < nAct);
	CONT(abpF.rep.facStf == 0. && !(nlPnt[0] < nAct && nlPnt[1] < nAct));
	mode = 0;
	for (k = 0; k < 2; k++) {
		if (nlPnt[k] < nAct) {
			ind[2 * k] = iAct[P2A(act.cyl.l,nlPnt[k],0,2)];
			ind[2 * k + 1] = iAct[P2A(act.cyl.l,nlPnt[k],1,2)];
		}
		else {
			ind[2 * k] = iAbp[nlPnt[k] - nAct];
			ind[2 * k + 1] = ind[2 * k];
			mode++;
		}
	}
	CONT(ind[0] < 0 || ind[1] < 0 || ind[2] < 0 || ind[3] < 0);
	if (mode == 0) {
		CONT(ind[0] == ind[2] || ind[1] == ind[3] 
				|| ind[0] == ind[3] || ind[1] == ind[2]);
		V2SET_ALL(dist, actF.dia);
	}
	else if (mode == 1) {
		V2SET(dist, actF.dia, abpF.repDia[K_ABP(ind[2])]);
	}
	else {
		if (ISMTF(K_ABP(ind[0]))) {
			CONT(P2A(abp.ch,ind[0],3,nChAb) == nlPnt[1] - nAct	
					|| P2A(abp.ch,ind[0],4,nChAb) == nlPnt[1] - nAct);
		}
		V2SET(dist, abpF.repDia[K_ABP(ind[0])], abpF.repDia[K_ABP(ind[2])]);
	}
	for (k = 0; k < 4; k++) {
		if (nlPnt[(int)(k / 2)] < nAct) {
			r_p[k] = &P2(act.r,ind[k],0);
			f_p[k] = &P2(act.f,ind[k],0);
		}
		else {
			r_p[k] = &P2(abp.r,ind[k],0);
			f_p[k] = &P2(abp.f,ind[k],0);
		}
	}
	


	if (mode == 1) {
		if (ISMTF(K_ABP(ind[2]))) {
			abpInd = P2A(abp.ch,ind[2],3,nChAb);
			if (abpInd > -1) {
				if (iAbp[abpInd] > -1) {
					r_p[3] = &P2(abp.r,iAbp[abpInd],0);
					f_p[3] = &P2(abp.f,iAbp[abpInd],0);
				
					for(k = 0; k < 4; k++) {
						V3COPY(rPnt[k], r_p[k]);
					}
					dist2 = AVG2(dist[0], dist[1]);
					CalcRepulsiveForcesSubSubSubroutine(rPnt, 1, 3);
					len = CalcSegSegDist(&rPnt[0], &rPnt[2], dr, ratio, 0);
					CONT(len == 0.);
					if (len < dist2) {	
						f = abpF.rep.stf * (dist2 - len);
						f /= len;
						VS3COPY(fi, dr, f);
						VVS3ADD(f_p[0], fi, REVSIGN(1. - ratio[0]));
						VVS3ADD(f_p[1], fi, REVSIGN(ratio[0]));
						VVS3ADD(f_p[2], fi, 1. - ratio[1]);
						VVS3ADD(f_p[3], fi, ratio[1]);
					}
				}
			}
		}
		else {
			for(k = 0; k < 4; k++) {
				V3COPY(rPnt[k], r_p[k]);
			}
			dist2 = AVG2(dist[0], dist[1]);
			CalcRepulsiveForcesSubSubSubroutine(rPnt, 1, 2);
			len = CalcSegPntDist(&rPnt[0], rPnt[2], dr, ratio);
			CONT(len == 0.);
			if (len < dist2) {	
				f = abpF.rep.stf * (dist2 - len);
				f /= len;
				VS3COPY(fi, dr, f);
				VVS3ADD(f_p[0], fi, REVSIGN(1. - ratio[0]));
				VVS3ADD(f_p[1], fi, REVSIGN(ratio[0]));
				VV3ADD(f_p[2], fi);
			}
		}
	}	
	else if (mode == 2) {
		if ((!ISMTF(K_ABP(ind[2])) && ISMTF(K_ABP(ind[0]))) 
				|| (ISMTF(K_ABP(ind[2])) && !ISMTF(K_ABP(ind[0])))) {
			if (ISMTF(K_ABP(ind[0]))) {
				ind2[0] = 0;
				ind2[1] = 2;
			}
			else {
				ind2[0] = 2;
				ind2[1] = 0;
			}
			abpInd = P2A(abp.ch,ind[ind2[0]],3,nChAb);
			CONT(abpInd < 0);
			CONT(iAbp[abpInd] < 0);
			r_p[ind2[0] + 1] = &P2(abp.r,iAbp[abpInd],0);
			f_p[ind2[0] + 1] = &P2(abp.f,iAbp[abpInd],0);
					
				for(k = 0; k < 4; k++) {
				V3COPY(rPnt[k], r_p[k]);
			}
			dist2 = AVG2(dist[0], dist[1]);
			CalcRepulsiveForcesSubSubSubroutine(rPnt, 1, 3);
			len = CalcSegPntDist(&rPnt[ind2[0]], rPnt[ind2[1]], dr, ratio);
			CONT(len == 0.);
			if (len < dist2) {	
				f = abpF.rep.stf * (dist2 - len);
				f /= len;
				VS3COPY(fi, dr, f);
				VVS3ADD(f_p[ind2[0]], fi, REVSIGN(1. - ratio[0]));
				VVS3ADD(f_p[ind2[0] + 1], fi, REVSIGN(ratio[0]));
				VV3ADD(f_p[ind2[1]], fi);
			}
		}
		else if (ISMTF(K_ABP(ind[2])) && ISMTF(K_ABP(ind[0]))) {
			abpInd = P2A(abp.ch,ind[2],3,nChAb);
			abpInd2 = P2A(abp.ch,ind[0],3,nChAb);
			if (abpInd > -1 && abpInd2 > -1) {
				if (iAbp[abpInd] > -1 && iAbp[abpInd2] > -1) {
					r_p[1] = &P2(abp.r,iAbp[abpInd2],0);
					f_p[1] = &P2(abp.f,iAbp[abpInd2],0);
					r_p[3] = &P2(abp.r,iAbp[abpInd],0);
					f_p[3] = &P2(abp.f,iAbp[abpInd],0);
					
					for(k = 0; k < 4; k++) {
						V3COPY(rPnt[k], r_p[k]);
					}
					dist2 = AVG2(dist[0], dist[1]);
					CalcRepulsiveForcesSubSubSubroutine(rPnt, 1, 3);
					len = CalcSegSegDist(&rPnt[0], &rPnt[2], dr, ratio, 0);
					CONT(len == 0.);
					if (len < dist2) {	
						f = abpF.rep.stf * (dist2 - len);
						f /= len;
						VS3COPY(fi, dr, f);
						VVS3ADD(f_p[0], fi, REVSIGN(1. - ratio[0]));
						VVS3ADD(f_p[1], fi, REVSIGN(ratio[0]));
						VVS3ADD(f_p[2], fi, 1. - ratio[1]);
						VVS3ADD(f_p[3], fi, ratio[1]);
					}
				}
			}
		}
		else {
			for(k = 0; k < 4; k++) {
				V3COPY(rPnt[k], r_p[k]);
			}
			dist2 = AVG2(dist[0], dist[1]);
			CalcRepulsiveForcesSubSubSubroutine(rPnt, 2, 2);
			len = CalcVecDist(dr, rPnt[0], rPnt[2], 0);
			CONT(len == 0.);
			if (len < dist2) {	
				f = abpF.rep.stf * (dist2 - len);
				f /= len;
				VS3COPY(fi, dr, f);
				VV3SUB(f_p[0], fi);
				VV3ADD(f_p[2], fi);
			}
		}
	}
	else {
		// rPnt[][] has the information about positions of both ends of two 
		// cylindrical segments. If two end points of one cylindrical segment
		// in r_p[] are connected through a periodic boundary, rPnt[] has
		// offset position to calculate a minimum distance between two segments.
	
		// Adjust segment positions
		for(k = 0; k < 4; k++) {
			V3COPY(rPnt[k], r_p[k]);
		}
		CalcRepulsiveForcesSubroutine(ind, nlPnt, rPnt, f_p, dist, mode);	
	}
  }
}

/*------------------------ Related to repulsive forces -----------------------*/

/*------------------------- Related to spring forces -------------------------*/

// Subroutine for CalcSpringForces().
// In a case of motors, two spring forces are applied in directions 
// longitudinal and perpendicular directions to the axis of actin filaments.
void CalcSpringForcesSubroutine(int ind1, int ind2, int loc, 
		int *chkBondForce) {
  double dr[NDIM], fi[NDIM], f, f2, len, lenR, ratio, rPos[NDIM], drAxis[NDIM];
  double rPos2[NDIM], dr2[NDIM], len2, fac, fAll;
  int k, CS, side, kind;
  int locInd1, locInd2, locNextActInd, chkErr;

  CS = 1;
  locInd1 = iAct[ind1];
  // Avoid double-counting bonds between actin segments
  if (loc < 2 && locInd1 < nActMe) {
	if (P2A(chkBondForce,locInd1,loc,2) != 0) { CS = 0; }
  }
  if (CS == 1) {
	// actin-actin bond
	if (loc < 2) { 
		locInd2 = iAct[ind2];
		len = CalcVecDist(dr, &P2(act.r,locInd1,0), &P2(act.r,locInd2,0), 0); 
        if (len * 1.4 > maxActCh) {
			maxActCh = 1.4 * len; 
		}
		f = SPRING(actF.spr.stf, len, actF.spr.eq); 

		if (confVmdInfo.tgl != 0) {
			if (locInd1 < nActMe) { 
				recAct.sprF[locInd1] += REVSIGN(f); 
			}
			if (locInd2 < nActMe) {
				recAct.sprF[locInd2] += REVSIGN(f);
			}
		}
		// Check possible erros
		chkErr = CheckLargeForce(f, 4);
		if (chkErr != 1) {
			RecordErrorSpringForce(ind1, ind2, f, len, 0);
		}
		AddSpringForce(f, len, dr, &P2(act.f,locInd1,0), &P2(act.f,locInd2,0));

		if (locInd2 < nActMe) { 
			P2A(chkBondForce,locInd2,1 - loc,2) = 1; 
		}
	}
	// actin-ABP bond
	else { 
		locNextActInd = iAct[P2A(act.ch,locInd1,0,nChAc)];
		locInd2 = iAbp[ind2];
		kind = K_ABP(locInd2);
		// Find the location of ABP on the actin segment.
		ratio = (double)((int)((loc - 2) / nChAcY)) / (double)nChAcX;
		CalcPosOnActSeg(&P2(act.r,locInd1,0), &P2(act.r,locNextActInd,0), 
				rPos, ratio);
		len = CalcVecDist(dr, rPos, &P2(abp.r,locInd2,0), 0); 
		if (ISMTF(kind)) {		
			CalcVec(drAxis, &P2(act.r,locInd1,0), &P2(act.r,locNextActInd,0));
			fac = V3DOT(dr, drAxis) / V3LEN_SQ(drAxis);
			VS3SUB(rPos2, rPos, drAxis, fac);
			ApplyBoundCondVector(rPos2, -1, 0);
	 		len = CalcVecDist(dr, rPos2, &P2(abp.r,locInd2,0), 0);
			len2 = CalcVecDist(dr2, rPos, rPos2, 0);
		}
		side = (P2A(abp.ch,locInd2,0,nChAb) == ind1) ? 0 : 1;
		f = SPRING(abpF.spr[kind].stf, len, abpF.spr[kind].eq);
		if (ISMTF(kind)) {
			f2 = SPRING(abpF.spr[kind].stf2, len2, 0.);
		}
		fAll = (ISMTF(kind)) ? sqrt(SQR(f) + SQR(f2)) : f;
		// Record force for coloring method via VMD
		if (confVmdInfo.tgl != 0) {
			if (locInd1 < nActMe) { 
				recAct.sprF[locInd1] += REVSIGN(fAll * (1 - ratio)); 
			}
			if (ratio > 0) {
				if (locNextActInd < nActMe) {	
					recAct.sprF[locNextActInd] += REVSIGN(fAll * ratio); 
				}
			}
			if (locInd2 < nAbpMe) {
				recAbp.sprF[locInd2] += REVSIGN(fAll);
			}
		}
		// Check possible errors
		chkErr = CheckLargeForce(fAll, 5);
		if (chkErr != 1) {
			RecordErrorSpringForce(ind1, ind2, fAll, len, 1);
		}
		// Distribute the calculated forces
		f /= len;
		if (ISMTF(kind) && len2 > 0) { f2 /= len2; }
		FOR_NDIM(k) {
			fi[k] = f * dr[k];
			if (ISMTF(kind)) { fi[k] += f2 * dr2[k]; }
	  		P2(act.f,locInd1,k) += fi[k] * (1 - ratio);
 			if (ratio > 0) {
	  			P2(act.f,locNextActInd,k) += fi[k] * ratio;
			}
			P2(abp.f,locInd2,k) -= fi[k];
		}
		// A positive number is tension, and a negative number is compression.
		P2A(recInstSprFabp,locInd2,side * 2,4) = REVSIGN(f * len);
  		if (recLongF.tgl != 0 || nMot > 0) {
			// Calculate the logitudinal component of a force acting on the 
			// motor arm toward pointed ends.
			if (ISMTF(kind)) {
				P2A(recInstSprFabp,locInd2,side * 2 + 1,4) 
							= f2 * len2 * ((fac < 0) ? -1. : 1.);
			}
			else {
				if (f < 0) {
					CalcUnitVec(drAxis, &P2(act.r,locInd1,0), 
							&P2(act.r,locNextActInd,0));
					P2A(recInstSprFabp,locInd2,side * 2 + 1,4) 
							= V3DOT(fi, drAxis);
				}
				else {
					P2A(recInstSprFabp,locInd2,side * 2 + 1,4) = 0;
				}
			}
		}
	}
  }
}

// This function scans the array "act.ch" and "abp.ch" to look for bonds between
// actin and actin or between actin and ABP. If found, it computes the spring 
// force. "chkActBF" precludes the double-counting of actin-actin bonds.
void CalcSpringForces(void) {
  int n, k, ind1, ind2, ind3, side;
  int *chkActBF;
  MALLOC(chkActBF,int,nActMe*2);
  memset(chkActBF, 0, sizeof(int) * nActMe * 2);

  FOR_ACTME(n) {
	CONT(ISACTM(n));
	for(k = 0; k < nChAc; k++) {
		ind2 = P2A(act.ch,n,k,nChAc);
		CONT(ind2 < 0);
 		CalcSpringForcesSubroutine(act.id[n], ind2, k, chkActBF);
	}
  }
  // If actin belongs to adjacent subdomain, and ABP belongs to the current
  // subdomain, it was neglected in the above calculation, so it has to done
  // here.
  FOR_ABPME(n) {
	for(k = 0; k < 2; k++) {
		ind1 = P2A(abp.ch,n,k,nChAb);
		CONT(ind1 < 0);
		CONT(iAct[ind1] < nActMe);
		side = FindAbpActinChain(iAct[ind1], abp.id[n], 0);
		CalcSpringForcesSubroutine(ind1, abp.id[n], side, 
				chkActBF);
	}
  }

  if (nChAcX > 1) {
	FOR_ACTCP(n) {
		ind1 = n + nActMe;
	    ind2 = P2A(act.ch,ind1,0,nChAc);
		CONT(ind2 < 0);
		CONT(iAct[ind2] < 0 || iAct[ind2] >= nActMe);
		for(k = 2; k < nChAc; k++) {
		    ind3 = P2A(act.ch,ind1,k,nChAc);
			CONT(ind3 < 0);
			CONT(iAbp[ind3] < nAbpMe);
			CalcSpringForcesSubroutine(act.id[ind1], ind3, k, chkActBF);
		}
	}
  }	
  free(chkActBF);
}

/*------------------------- Related to spring forces -------------------------*/

/*----------------------------- Related to actins ----------------------------*/

// Subroutine for CalcFilaBendForces().
void CalcFilaBendForcesSubroutine(int ind1, int ind2, int ind3) {
  double dr1[NDIM], dr2[NDIM], f, f1[NDIM], f2[NDIM], fSum[NDIM];
  int chkErr;

  CalcVec(dr1, &P2(act.r,ind1,0), &P2(act.r,ind2,0));
  CalcVec(dr2, &P2(act.r,ind3,0), &P2(act.r,ind1,0));	

  f = CalcBendForce(dr1, dr2, &P2(act.f,ind2,0), &P2(act.f,ind1,0),
		&P2(act.f,ind3,0), actF.bend.stf, actF.bend.eq, f1, f2);
  if (confVmdInfo.tgl != 0) {
	V3ADD(fSum, f1, f2);
	if (ind2 < nActMe) { recAct.bendF[ind2] += V3LEN(f1); }
	if (ind1 < nActMe) { recAct.bendF[ind1] += V3LEN(fSum); }
	if (ind3 < nActMe) { recAct.bendF[ind3] += V3LEN(f2); }
  }
  // Check possible errors
  chkErr = CheckLargeForce(f, 6);
  if (chkErr != 1) { 
	RecordErrorBendingForce(act.id[ind1], act.id[ind2],
			act.id[ind3], f, 0);
  }
}
// Calculate bending forces keeping actin filaments straight.
void CalcFilaBendForces(void) {
  int n, *pArr;

  FOR_ACTME(n) {
	pArr = &P2A(act.ch,n,0,nChAc);
	CONT(!(pArr[0] > -1 && pArr[1] > -1));
	CalcFilaBendForcesSubroutine(n, iAct[pArr[0]], iAct[pArr[1]]);
  }
  // If two of the three actin segments belong to adjacent subdomain, and
  // if the other belongs to the current subdomain, it is neglected in the
  // above calculation. So, it has to be performed here.
  for(n = 0; n < nActCp; n++) {
	pArr = &P2A(act.ch,n + nActMe,0,nChAc);
	CONT(!(pArr[0] > -1 && pArr[1] > -1));
	CONT(!((iAct[pArr[0]] < nActMe || iAct[pArr[1]] < nActMe) 
			&& (iAct[pArr[0]] > -1 && iAct[pArr[1]] > -1)));
	CalcFilaBendForcesSubroutine(n + nActMe, iAct[pArr[0]], iAct[pArr[1]]);
  }
}

/*----------------------------- Related to actins ----------------------------*/

/*------------------------------ Related to ABPs -----------------------------*/

void CalcFilaBranchForces(void) {
  int n, k, *pArr, loc;
  int actInd[2], locActInd[2], nextActInd[2], locNextActInd[2];
  double dr[5][NDIM], ang, f, rPos[NDIM], crs[3][NDIM], mag[2], ff;

  actBch.stf = KB_NM2S(1e-18);

  FOR_ABPME(n) {
	pArr = &P2A(abp.ch,n,0,nChAb);
	CONT(!(pArr[0] > -1 && pArr[1] > -1 && K_ABP(n) == 0));
	for(k = 0; k < 2; k++) {
		actInd[k] = pArr[k];
		locActInd[k] = iAct[actInd[k]];
		nextActInd[k] = P2A(act.ch,locActInd[k],0,nChAc);
		locNextActInd[k] = iAct[nextActInd[k]];
		mag[k] = CalcVecDist(dr[k], &P2(act.r,locNextActInd[k],0), 
				&P2(act.r,locActInd[k],0), 0);
	}
	CONT(!(mag[0] > POS_SMALL_VALUE && mag[1] > POS_SMALL_VALUE));

	ang = V3ANG(dr[0], dr[1]);
	f = actBch.stf * (ang - DEG2RAD(angFila));
	if (ang < POS_SMALL_VALUE) {
		CalcPosOnActSegSide(&P2(act.r,locActInd[0],0), 
				&P2(act.r,locNextActInd[0],0), rPos, loc);
		CalcVec(dr[2], &P2(act.r,locActInd[1],0), rPos);
	}
	else {
		V3CROSS(dr[2], dr[0], dr[1]);
	}
	for(k = 0; k < 2; k++) {
		ff = f / mag[k];
		V3CROSS(dr[k + 3], dr[k], dr[2]);
		NormVec(dr[k + 3]);
		VVS3ADD(&P2(act.f,locActInd[k],0), dr[k + 3], (k == 0 ? 1. : -1.) * ff);
		VVS3ADD(&P2(act.f,locNextActInd[k],0), dr[k + 3], 
				(k == 0 ? -1. : 1.) * ff);
	}
	loc = FindAbpActinChain(locActInd[0], abp.id[n], 0);	
	CalcPosOnActSegSide(&P2(act.r,locActInd[0],0), &P2(act.r,locNextActInd[0],0), 
			rPos,  loc);
	CalcVec(dr[2], &P2(act.r,locActInd[1],0), rPos);
	V3CROSS(crs[0], dr[2], dr[0]);
	V3CROSS(crs[1], dr[2], dr[1]);
	ang = V3ANG(crs[0], crs[1]);
	V3CROSS(crs[2], crs[0], crs[1]);
	if (V3DOT(crs[2], dr[2]) < 0.) {
		ang *= -1.;
	}	
	f = actBch.stf * ang;
	NormVec(crs[0]);
	NormVec(crs[1]);
	ff = f / mag[0];
	VVS3ADD(&P2(act.f,locNextActInd[0],0), crs[0], ff);
	VVS3ADD(&P2(act.f,locActInd[0],0), crs[0], -1. * ff);
	ff = f / mag[1];
	VVS3ADD(&P2(act.f,locNextActInd[1],0), crs[1], -1. * ff);
	VVS3ADD(&P2(act.f,locActInd[1],0), crs[1], ff);
  } 
}

// 1st subroutine for CalcAbpBendForces(), which calculates bending force
// maintaining right angle between filament axis and the arm of ABP
void CalcAbpBendForcesSubroutine1(int abpInd, int actInd) {
  int k, locActInd, locAbpInd, nextActInd, locNextActInd, chkErr;
  double f, fi[4][NDIM], fAdd[3][NDIM], ratio, fSum;
  double dr1[NDIM], dr2[NDIM], dr3[NDIM];
  
  locActInd = iAct[actInd];
  locAbpInd = iAbp[abpInd];
  nextActInd = P2A(act.ch,locActInd,0,nChAc);
  locNextActInd = iAct[nextActInd];

  if (K_ABP(locAbpInd) == 0 && P2A(abp.ch,locAbpInd,1,nChAb) == actInd) {
	CalcVec(dr1, &P2(act.r,locActInd,0), &P2(abp.r,locAbpInd,0));
	CalcVec(dr2, &P2(act.r,locNextActInd,0), &P2(act.r,locActInd,0));	

	f = CalcBendForce(dr1, dr2, &P2(abp.f,locAbpInd,0), &P2(act.f,locActInd,0),
			&P2(act.f,locNextActInd,0), abpF.a90[0].stf, DEG2RAD(90. - angFila), 
			fAdd[0], fAdd[1]);
	if (confVmdInfo.tgl != 0) {
		V3ADD(fAdd[2], fAdd[0], fAdd[1]);
		if (locAbpInd < nAbpMe) { recAbp.bendF[locAbpInd] += V3LEN(fAdd[0]); }
		if (locActInd < nActMe) { recAct.bendF[locActInd] += V3LEN(fAdd[2]); }
		if (locNextActInd < nActMe) { recAct.bendF[locNextActInd] 
				+= V3LEN(fAdd[1]); }
	}
	// Check possible errors
	chkErr = CheckLargeForce(f, 8);
	if (chkErr != 1) { 
		RecordErrorBendingForce(abpInd, actInd, nextActInd, f, 3); 
	}
  }
  else {
	ratio = CalcVecActinAbp(dr1, actInd, abpInd, 0);
	CalcVec(dr2, &P2(act.r,locNextActInd,0), &P2(act.r,locActInd,0));

	for(k = 0; k < 2; k++) {
		if (k == 0) { VS3COPY(dr3, dr2, 1. - ratio); }
		else { 
			if (ratio == 0.) {
				V3SET_ALL(fi[2], 0.);
				V3SET_ALL(fi[3], 0.);
				continue;
			}
			VS3COPY(dr3, dr2, REVSIGN(ratio)); 
		}
		f = CalcBendForceSubroutine(dr1, dr3, abpF.a90[K_ABP(locAbpInd)].stf, 
				abpF.a90[K_ABP(locAbpInd)].eq, fi[k * 2], fi[k * 2 + 1]);
	}
	FOR_NDIM(k) {
		fAdd[0][k] = fi[0][k] + fi[2][k];
		fSum = fi[0][k] + fi[1][k] + fi[2][k] + fi[3][k];
		fAdd[1][k] = fi[3][k] - fSum * (1. - ratio);
		fAdd[2][k] = fi[1][k] - fSum * ratio;
	}
	VV3ADD(&P2(abp.f,locAbpInd,0), fAdd[0]);
	VV3ADD(&P2(act.f,locActInd,0), fAdd[1]);
	VV3ADD(&P2(act.f,locNextActInd,0), fAdd[2]);
	// Record forces for coloring method via VMD
	if (confVmdInfo.tgl != 0) {
		if (locAbpInd < nAbpMe) { recAbp.bendF[locAbpInd] += V3LEN(fAdd[0]); }
		if (locActInd < nActMe) { recAct.bendF[locActInd] += V3LEN(fAdd[1]); }
		if (locNextActInd < nActMe) { 
			recAct.bendF[locNextActInd] += V3LEN(fAdd[2]); 
		}
	}
	// Check errors
	chkErr = CheckLargeForce(f, 8);
	if (chkErr != 1) { 
		RecordErrorBendingForce(abpInd, actInd, nextActInd, f, 3); 
	}
  }
}

// 2nd subroutine for CalcAbpBendForces(), which calculates bending forces
// for angle formed by actin-ABP-actin.
void CalcAbpBendForcesSubroutine2(int abpInd, int *actInd) {
  int m, k, chkErr, locAbpInd, locActInd[2], locNextActInd[2];
  double mag, f, fSum[NDIM], f2[2][NDIM], dr[2][NDIM], ratio[2];

  locAbpInd = iAbp[abpInd];
  for(k = 0; k < 2; k++) {
	locActInd[k] = iAct[actInd[k]];
	ratio[k] = CalcVecActinAbp(dr[k], actInd[k], abpInd, 0);
	locNextActInd[k] = iAct[P2A(act.ch,locActInd[k],0,nChAc)];
  }
  V3REVSIGN(dr[0]);

  f = CalcBendForceSubroutine(dr[0], dr[1], abpF.bend[K_ABP(locAbpInd)].stf, 
		abpF.bend[K_ABP(locAbpInd)].eq, f2[0], f2[1]);
  
  FOR_NDIM(k) {
	P2(abp.f,locAbpInd,k) -= f2[0][k] + f2[1][k];
	for(m = 0; m < 2; m++) {
		P2(act.f,locActInd[m],k) += f2[m][k] * (1. - ratio[m]);
		P2(act.f,locNextActInd[m],k) += f2[m][k] * ratio[m];
	}
  }
  if (confVmdInfo.tgl != 0) {
	V3ADD(fSum, f2[0], f2[1]);
	for(k = 0; k < 2; k++) {
		mag = V3LEN(f2[k]);
		if (locActInd[k] < nActMe) { 
			recAct.bendF[locActInd[k]] += mag * (1. - ratio[k]); 
		}
		if (locNextActInd[k] < nActMe) { 
			recAct.bendF[locNextActInd[k]] += mag * ratio[k]; 
		}
	}
	if (locAbpInd < nAbpMe) { recAbp.bendF[locAbpInd] += V3LEN(fSum); }
  }
  chkErr = CheckLargeForce(f, 7);
  if (chkErr != 1) { 
	RecordErrorBendingForce(abpInd, actInd[0], actInd[1], f, 1);
  }
}

// Calculate two kinds of bending forces related to ABP
void CalcAbpBendForces(void) {
  int n, k, kind, abpInd, locAbpInd;
  int *pArr, locActInd[2], locNextActInd[2], CS[2];

  FOR_ABPME(n) {
	CONT(ISABPM(n));
	pArr = &P2A(abp.ch,n,0,nChAb);
	kind = K_ABP(n);
	CONT(ISMTF(kind));
	if (abpF.a90[kind].stf > 0.) {
		for(k = 0; k < 2; k++) {
			CONT(pArr[k] < 0);
		    CalcAbpBendForcesSubroutine1(abp.id[n], pArr[k]);
		}
	}
	if (abpF.bend[kind].stf > 0.) {
		if (pArr[0] > -1 && pArr[1] > -1) {
			CalcAbpBendForcesSubroutine2(abp.id[n], &pArr[0]);
		}
	}
  }

  // If ABP belongs to the adjacent subdomain, but if actin which is involved
  // with calculation of these bending forces belongs to the current subdomain,
  // the calculation has to be done separately. 
  for(n = 0; n < nAbpCp; n++) {
	locAbpInd = n + nAbpMe;
	abpInd = abp.id[locAbpInd];
	pArr = &P2A(abp.ch,locAbpInd,0,nChAb);
	CONT(ISMTF(pArr[2]));

	V2SET_ALL(CS, 0);
	for(k = 0; k < 2; k++) {
		CONT(pArr[k] < 0);
		locActInd[k] = iAct[pArr[k]];
		CONT(locActInd[k] < 0);
		locNextActInd[k] = iAct[P2A(act.ch,locActInd[k],0,nChAc)];
		CONT(locNextActInd[k] < 0);
		CONT(!(locActInd[k] < nActMe || locNextActInd[k] < nActMe));
		CS[k] = 1;
	}
	if (abpF.a90[pArr[2]].stf > 0.) {
		for(k = 0; k < 2; k++) {
			CONT(CS[k] != 1);
			CalcAbpBendForcesSubroutine1(abpInd, pArr[k]);
		}
	}
	if (abpF.bend[pArr[2]].stf > 0.) {
		CONT(!(CS[0] == 1 && CS[1] == 1));
		CalcAbpBendForcesSubroutine2(abpInd, pArr);
	}
  }
}

void CalcMotorBackboneForces(void) {
}

/*------------------------------ Related to ABPs -----------------------------*/

/*-------------------------- Tools for calculation ---------------------------*/

// Calculate cosine from two vectors.
void CalcCosine(double *d1, double *d2, double *cc11, double *cc12, 
		double *cc22, double *ccD, double *cc) {
  *cc11 = V3LEN_SQ(d1);
  *cc12 = V3DOT(d1, d2);
  *cc22 = V3LEN_SQ(d2);
  if (*cc11 < POS_SMALL_VALUE || *cc22 < POS_SMALL_VALUE) { return; }
  *ccD = sqrt (*cc11 * (*cc22));
  *cc = (*cc12) / (*ccD); // cosine
}

double CalcBendForceSubroutine(double *dr1, double *dr2, double stiff,
        double ang, double *f1, double *f2) {
  int k;
  double fe1[NDIM], fe2[NDIM], mag1, mag2, f, ff1, ff2;
  double c11, c12, c22, cD, c;

  CalcCosine(dr1, dr2, &c11, &c12, &c22, &cD, &c);

  f = stiff * (Acos(c) - ang);
 
  VS3SUB(fe2, dr1, dr2, c12 / c22);
  VS3SUB(fe1, dr2, dr1, c12 / c11);
  V3REVSIGN(fe1);
  mag1 = sqrt(c11 * V3LEN_SQ(fe1));
  mag2 = sqrt(c22 * V3LEN_SQ(fe2));
  if (mag1 > POS_SMALL_VALUE && mag2 > POS_SMALL_VALUE) {
	ff1 = f / mag1;
	ff2 = f / mag2;
	FOR_NDIM(k) {
		f1[k] = ff1 * fe1[k];
		f2[k] = ff2 * fe2[k];
	}
  }
  else {
    V3SET_ALL(f1, 0.);
    V3SET_ALL(f2, 0.);
  }
  return f;
}

// Subroutine for other subroutines calculating bending forces.
double CalcBendForce(double *dr1, double *dr2, double *force1, double *force2, 
		double *force3, double stiff, double ang, double *f1, double *f2) {
  int k;
  double f;

  f = CalcBendForceSubroutine(dr1, dr2, stiff, ang, f1, f2);

  FOR_NDIM(k) {
	force1[k] += f1[k];
	force2[k] -= f1[k] + f2[k];
	force3[k] += f2[k];
  }
  return f;
}

void AddSpringForce(double f, double len, double *dr, double *f1, double *f2) {
  double fi[NDIM];

  f /= len;
  VS3COPY(fi, dr, f);
  VV3ADD(f1, fi);
  VV3SUB(f2, fi);
}

/*-------------------------- Tools for calculation ---------------------------*/
