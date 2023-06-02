/* estimNe.c - NEUTRAL MODEL */

/* ********************************************************************* */

#include "libhdr"
#include "ranlib.h"
#define NN 5000   /* max number of NIND */
#define GG 1000   /* max number of gen */
#define MM 1001   /* max number of NCRO */
#define maxmpt 5001
#define normalthreshold 30  

int NIND, NINDs, NINDNP, bott, bottBP, bottL, bottr, botte, bottLe, ib, bt, rate, maxexp;
int init, tLINES, genBPr, genBPb, endbm, genBPe, endbe, genLe, endLe, endmL, expns[NN], start, end; 
int CRO, CROX, KCRO, NCRO, NLOCI, TOTLOCI, NSEGLOCNP;
int gen, tgen, i, ran_i, j, k, ran_k, ran_h, l, ran_l, m, rep, replicates;
int a, b, expr, expBP, subexp, nsubexp, lin, lines, getfiles[GG];
int muts, r, mutants_ocurred, lastinmutantspoissontable, lastinrecombinantpoissontable;
int n, ss, countNoSS, SNP[MM][31], RM[NN], NoSS_k[60001], NoSS_l[60001];
int father[NN], fatherBP[NN], fatherBPL[NN], fatherL1[NN], fatherL2[NN], fatherEXP1[NN];
int mother[NN], motherBP[NN], motherBPL[NN], motherL1[NN], motherL2[NN], motherEXP1[NN];
int g[NN][MM][2], gNP[10001][MM][2], gBP[NN][MM][2], gBPL[NN][MM][2], gdf[NN][MM][2], sp[NN][MM][2];
int gL1[NN][MM][2], gL2[NN][MM][2], gEXP1[NN][MM][2];
int initialgen[MM][31], initialgenNP[MM][31], initialgenBP[MM][31], initialgenBPL[MM][31];
int initialgenL1[MM][31], initialgenL2[MM][31], initialgenEXP1[MM][31];
int NOFF, bottOFF, spmother[NN], spfather[NN], matingsystem, isolation;

double experiment, migration, mig, L, Lambda, prop, AA, Aa, aa, sk2, ne;
double q[MM][31], qNP[MM][31], qBP[MM][31], qBPL[MM][31], qL1[MM][31], qL2[MM][31], qEXP1[MM][31], meanQlines[GG][MM][31];
double mat[NN][NN], smat[NN][NN], smatBP[NN][NN], smatBPL[NN][NN], smatL1[NN][NN], smatL2[NN][NN], smatEXP1[NN][NN];
double mutantspoissontable[maxmpt], recombinantpoissontable[maxmpt];
double H, Hw[8][GG], Hw_btmean[8][GG], HO_btmean[8][GG], sumq[8][GG], sumq2[8][GG], Fstjunk;
double Qmean[8][GG], Qbottle[8][GG], HTjk[8][GG], HTjunklines[GG], Vq_btmean[8][GG], Vq[8][GG]; 

struct acc N[8][GG], nbott[8][GG], SK2[8][GG], fp[8][GG], Fp[8][GG], Fst[8][GG];
struct acc q_00[8][101], q_00_01[8][101], q_01_02[8][101], q_02_03[8][101], q_03_04[8][101], q_04_05[8][101];
struct acc q_05_06[8][101], q_06_07[8][101], q_07_08[8][101], q_08_09[8][101], q_09_10[8][101], q_10[8][101];
struct acc Hwglobal[8][GG], HOglobal[8][GG], Hwbottle[8][GG], HObottle[8][GG], HT[8][GG], HTlines[GG];
struct acc Qglobal[8][GG], meanQbottle[8][GG], lociseg[8][GG], Varq[8][GG], Varq_bt[8][GG];
struct acc Ne[8][GG], NeH[8][GG], NeH_bottle[8][GG], NeFst[8][GG], Fst_lines[GG], NeFst_lines[GG], Nevq[8][GG], Nevq_bottle[8][GG]; 

FILE *fptr, *fgenBP, *fgenL, *fgenLines, *fgenmL, *fdat, *fpop, *fsum; 
FILE *fdataBP, *fdataL1, *fdataEXP1R, *fdataEXP1B, *fdataEXP1E, *fdataEXP2LE, *fdataEXP3ML;
FILE *fmapBP, *fmapL1, *fmapEXP1R, *fmapEXP1B, *fmapEXP1E, *fmapEXP2LE, *fmapEXP3ML;

/* ********************************************************************* */

main()
{
	getinputs();

	if ((experiment!=0.2)&&(experiment!=2)&&(experiment!=3)) 	fgenBP = fopen ("genfile_BP.dat","w");
	if ((experiment==0.2)||(experiment==2)||(experiment==99))	fgenL = fopen ("genfile_Lines.dat","w");
	if ((experiment==3)||(experiment==99))				fgenmL = fopen ("genfile_Cross.Line.dat","w");
	if (experiment==99)						fgenLines = fopen ("genfile_Fst_Lines.dat","w");

	if ((experiment==0.1)||(experiment==99))   { fdataBP = fopen ("BP.ped","w"); 		 fmapBP = fopen ("BP.map","w"); }
	if ((experiment==1.0)||(experiment==99))   { fdataEXP1R = fopen ("BPred.ped","w"); 	 fmapEXP1R = fopen ("BPred.map","w"); }
	if ((experiment==1.1)||(experiment==99))   { fdataEXP1E = fopen ("BPbot_exp.ped","w");	 fmapEXP1E = fopen ("BPbot_exp.map","w"); }
	if ((experiment==1.2)||(experiment==99))   { fdataEXP1B = fopen ("BPbot.ped","w");	 fmapEXP1B = fopen ("BPbot.map","w"); }
	if ((experiment==0.2)||(experiment==99))   { fdataL1 = fopen ("LINE.ped","w");	 	 fmapL1 = fopen ("LINE.map","w"); }
	if ((experiment==2)||(experiment==99))     { fdataEXP2LE = fopen ("LINEexp.ped","w");	 fmapEXP2LE = fopen ("LINEexp.map","w"); }
	if ((experiment==3)||(experiment==99))     { fdataEXP3ML = fopen ("CrossL.ped","w"); 	 fmapEXP3ML = fopen ("CrossL.map","w"); }

	fsum = fopen ("summary_outline.dat","w");
	if (tracelevel!=0) 	fptr = fopen ("dfilename.dat","w");

	sum_outline();	
	recombination_masks();
	natural_population();

	for (rep=1; rep<=replicates; rep++)
	{
		if (tracelevel!=0)	fprintf (fptr,"\n***********************\n***** Replicate %d *****\n***********************\n", rep);
		tgen = 10;

		for (lin=0; lin<=7; lin++)
		{
			//BASE POPULATION AND LINES
			if (lin==0)							BP();
			else if ((lin==1) && ((experiment==0.2)||(experiment>=2.0)))	{lines = 1; LINES();}
			else if ((lin==2) && (experiment>=3.0))				{lines = 2; LINES();}
			
			//EXPERIMENT A
			else if ((lin==3) && (((experiment>=1.0)&&(experiment<2.0))||(experiment==99)))	REDUCTION_N();
			else if ((lin==4) && ((experiment==1.2)||(experiment==99)))			BOTTLENECK();
			else if ((lin==5) && ((experiment==1.1)||(experiment==99)))			BOTandEXP();
			
			//EXPERIMENT B
			else if ((lin==6) && ((experiment==2.0)||(experiment==99)))	LINE_EXPANSION();

			//EXPERIMENT C
			else if ((lin==7) && ((experiment==3.0)||(experiment==99)))	MIXED_LINE();
			else goto labelLIN;

			for (gen=start; gen<end; gen++)
			{
				if (tracelevel!=0)	fprintf (fptr,"\n*** generation %d ***\n", gen);
				accum(&N[lin][gen], (double)NIND);
				accum(&nbott[lin][gen], (double)bott);

				poisson_tables();
//				mutation_neutral();
				neutral_genes();

//				if (tracelevel!=0)		dumpoffspringaftermutation();
				if (gen%tgen == 0)		distribution_q();

				coancestry_matrix();
				if (tracelevel!=0)		dumpparents();

				if (getfiles[gen]==1)
				{
					if (lin==0)		data_files(fmapBP, fdataBP);
					else if (lin==1)	data_files(fmapL1, fdataL1);
					else if (lin==3)	data_files(fmapEXP1R, fdataEXP1R);
					else if (lin==4)	data_files(fmapEXP1B, fdataEXP1B);
					else if (lin==5)	data_files(fmapEXP1E, fdataEXP1E);
					else if (lin==6)	data_files(fmapEXP2LE, fdataEXP2LE);
					else if (lin==7)	data_files(fmapEXP3ML, fdataEXP3ML);
				}	

				mating();
//				if (tracelevel!=0)		dumpoffspring();

				if ((lin==0)&&(gen==tLINES-1))
				{
					save(gBPL, motherBPL, fatherBPL, smatBPL, qBPL, initialgenBPL);
					if ((experiment==0.2)||(experiment==2)||(experiment==3))	goto labelLIN;
				}

				if ((lin==3)&&(experiment==1.0)&&(gen==genBPr-1))	goto labelEND;
			}

			if ((lin==0)&&(experiment==0.1))	goto labelEND;
			if ((lin==4)&&(experiment==1.2))	goto labelEND;
			if ((lin==5)&&(experiment==1.1))	goto labelEND;
			if ((lin==6)&&((experiment==2.0)))	goto labelEND;

			if (lin==0)		save(gBP, motherBP, fatherBP, smatBP, qBP, initialgenBP);
			else if (lin==1)	save(gL1, motherL1, fatherL1, smatL1, qL1, initialgenL1);
			else if (lin==2)	save(gL2, motherL2, fatherL2, smatL2, qL2, initialgenL2);
			else if (lin==3)	save(gEXP1, motherEXP1, fatherEXP1, smatEXP1, qEXP1, initialgenEXP1);

			labelLIN: /* end of line */;
		}

		labelEND: /* end of replicate */;
	}

	printout();
	distribution_out();
	writeseed();
}

/* ********************************************************************* */

getinputs()
{
	tracestart();
	getseed();
	getrealandskip("Experiment: mant(0.1; 0.2), A(1.0; 1.1; 1.2), B(2), C(3), ALL(99)",&experiment, 0.0, 99.0);
	getintandskip("Mating system (0 circular, 1 offspring pool, 2 complete panmixia):",&matingsystem,0,2);
	getrealandskip("Migration between bottles:",&migration,0.0,1.0);
	getintandskip("NIND per Bottle (max 80):",&ib,2,80);
	getintandskip("Number Bottles - BP (max 60):",&bottBP,1,60);
	getintandskip("Generations before experiments - BP:",&init,1,1001);
	getintandskip("Formation of lines:",&tLINES,0,init);
	getintandskip("Number Bottles - Lines (max 60):",&bottL,0,60);
	getintandskip("if experiment=0.1, maintain Pb bottles isolated after formation of lines (no 1, yes 0):",&isolation,0,1);
	getintandskip("(expA)   Bottles bottleneck - BP (max 60):",&bottr,0,60);
	getintandskip("(expA)   Iniciation recovery - BP (generation):",&genBPr,0,1000);
	getintandskip("(expA)   Stop recovery - BP (generation):",&genBPb,0,1000);
	getintandskip("(expA) Mantain size until generation:",&endbm,0,1000);
	getintandskip("(expA 1.1) Bottles expansion (max 60):",&botte,0,60);
	getintandskip("(expA 1.1) Expansion until generation:",&genBPe,0,1000);
	getintandskip("(expA 1.1) Mantain size until generation:",&endbe,0,1000);
	getintandskip("(expB)   Bottles expansion (max 60):",&bottLe,0,60);
	getintandskip("(expB)   Expansion until generation:",&genLe,0,1000);
	getintandskip("(expB)   Mantain size until generation:",&endLe,0,1000);
	getintandskip("(expC)   Mantain mixed line until generation:",&endmL,0,1000);
	getrealandskip("Length per chromosome in Morgans(99:FreeRecom) :",&L,0.0,99.0);
	getintandskip("Number of Chromosomes (min 1, max 50):",&CRO,1,50);
	getintandskip("Chromosome X (no 0, yes 1):",&CROX,0,1);
	getintandskip("NCRO (min 1, max 2000):",&NCRO,1,2000);
	getintandskip("NLOCI (min 2, max 30):",&NLOCI,2,30);
	TOTLOCI = NCRO * NLOCI;
	KCRO = NCRO / CRO;
	getrealandskip("Lambda (neutral model):",&Lambda,0.0,(double)infinity);
	getintandskip("Sample size for data files:",&NINDs,1,4800);
	getintandskip("Replicates:",&replicates,1,infinity);
}

/* ********************************************************************* */

recombination_masks ()
{
	for (l=0; l<NLOCI; l++)   RM[l]=pow(2.0,(double)l);
}

/* ********************************************************************* */

natural_population()
{
//COMMENT
//	if (tracelevel!=0)	fprintf(fptr,"\nCOMMENT - natural_population\n");

	int g0, g1, dinitialgen, dN;
	double dadds, dadda, ds, da, dhs, dha;

	/* ***** take genotypic values of natural population ***** */

	fpop=fopen("popfile","r");

	fscanf(fpop,"%d", &dN);
	NINDNP = dN;

	for (i=0; i<dN; i++)
	for (k=0; k<NCRO; k++)
	{
		fscanf(fpop,"%d%d", &g0, &g1);
		gNP[i][k][0] = g0;
		gNP[i][k][1] = g1;
	}

	fclose(fpop);

	/* **************** take effects of genes **************** */

	fdat=fopen("datafile","r");

	fscanf(fdat,"%lf%lf", &dadds, &dadda);

	for (k=0; k<NCRO; k++) 
	for (l=0; l<NLOCI; l++)
	{
		fscanf(fdat,"%lf%lf%lf%lf%d", &ds, &da, &dhs, &dha, &dinitialgen);
		
		initialgenNP[k][l] = dinitialgen;
	}

	/* ******** frequences in the natural population ********* */

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	{
		AA=0.0; Aa=0.0; aa=0.0;

		for (i=0; i<NINDNP; i++)
		{
			if (((gNP[i][k][0] & RM[l])==RM[l])&&((gNP[i][k][1] & RM[l])==RM[l]))	aa+=1.0;
	    		else if (((gNP[i][k][0] & RM[l])!=RM[l])&&((gNP[i][k][1] & RM[l])!=RM[l]))	AA+=1.0;
		     	else	Aa+=1.0;
		}

		qNP[k][l] = (aa/(double)NINDNP)+(Aa/(2.0*(double)NINDNP));
		if (qNP[k][l] > 0.0)	NSEGLOCNP ++;

//COMMENT
//		if (tracelevel!=0)	fprintf(fptr,"k=%d\tl=%d\tAA=%1.0f\tAa=%1.0f\taa=%1.0f\tqNP[k][l]=%f\tinitialgenNP[k][l]=%d\n",k,l,AA,Aa,aa,qNP[k][l],initialgenNP[k][l]);
	} 

	fclose(fdat);

//COMMENT
//		if (tracelevel!=0)	for (i=0; i<NINDNP; i++)	fprintf(fptr,"gNP[%d][0][0]=%d	gNP[%d][0][1]=%d\n", i, gNP[i][0][0], i, gNP[i][0][1]);

}

/* ********************************************************************* */

BP()
{
	if (tracelevel!=0)	fprintf (fptr,"\n******* BASE POPULATION *******\n");

	bott = bottBP;
	NIND = ib * bott;
	start = 0;
	end = init;

	for (gen=start; gen<end; gen++)
	{
		expns[gen] = 0;
		if ((gen==end-1)&&((experiment==0.1)||(experiment==99)))	getfiles[gen] = 1;
		else	getfiles[gen] = 0;
	}

	sample_NP();
}

/* ********************************************************************* */

LINES()
{
	if (tracelevel!=0)	fprintf (fptr,"\n******* LINE %d (FROM BP) *******\n", lines);

	bott = bottL;
	NIND = ib * bott;
	start = tLINES;
	end = init;

	for (gen=start; gen<end; gen++)
	{
		expns[gen] = 0;
		if ((gen==end-1)&&((experiment==0.2)||(experiment==99)))	getfiles[gen] = 1;
		else	getfiles[gen] = 0;
	}

	sample(bottBP, gBPL, motherBPL, fatherBPL, smatBPL, qBPL, initialgenBPL);
}

/* ********************************************************************* */

REDUCTION_N()
{
	if ((tracelevel!=0)&&(experiment==99))	fprintf (fptr,"\n******* EXPERIMENT A *******\n");
	if ((tracelevel!=0)&&(experiment==1.0))	fprintf (fptr,"\n******* EXPERIMENT A (REDUCTION N) *******\n");
	if ((tracelevel!=0)&&(experiment==1.1))	fprintf (fptr,"\n******* EXPERIMENT A (BOTTLENECK & EXPANSION) *******\n");
	if ((tracelevel!=0)&&(experiment==1.2))	fprintf (fptr,"\n******* EXPERIMENT A (BOTTLENECK) *******\n");

	bott = bottr;
	NIND = ib * bott;
	rate = (int)(bottBP-bott)/(genBPb-genBPr);
	maxexp = bottBP;
	start = init;
	end = genBPb;

	for (gen=start; gen<end; gen++)
	{
		if (gen<genBPr)	expns[gen] = 0;
		else		expns[gen] = 1;

		if ((gen==genBPr-1)&&((experiment==1.0)||(experiment==99))) 	getfiles[gen] = 1;
		else	getfiles[gen] = 0;
	}

	sample(bottBP, gBP, motherBP, fatherBP, smatBP, qBP, initialgenBP);
}

/* ********************************************************************* */

BOTTLENECK()
{
	if (tracelevel!=0)	fprintf (fptr,"\n** mantain population size after bottleneck **\n");

	start = genBPb;
	end = endbm + 1;

	for (gen=start; gen<end; gen++)
	{
		expns[gen] = 0;
		if (gen==end-1)	getfiles[gen] = 1;
		else		getfiles[gen] = 0;
	}

	takeindividuals(gEXP1, motherEXP1, fatherEXP1, smatEXP1, qEXP1, initialgenEXP1);
}

/* ********************************************************************* */

BOTandEXP()
{
	if (tracelevel!=0)	fprintf (fptr,"\n** expansion after bottleneck **\n");

	rate = (int)(botte-bott)/(genBPe-genBPb);
	maxexp = botte;
	start = genBPb;
	end = endbe + 1;

	for (gen=start; gen<end; gen++)
	{
		if (gen<genBPe)	expns[gen] = 1; 
		else		expns[gen] = 0;

		if (gen==end-1)	getfiles[gen] = 1;
		else		getfiles[gen] = 0;

	}

	takeindividuals(gEXP1, motherEXP1, fatherEXP1, smatEXP1, qEXP1, initialgenEXP1);
}

/* ********************************************************************* */

LINE_EXPANSION()
{
	if (tracelevel!=0)	fprintf (fptr,"\n******* EXPERIMENT B (LINE EXPANSION) *******\n");

	bott = bottL;
	NIND = ib * bott;
	rate = (int)(bottLe-bott)/(genLe-init);
	maxexp = bottLe;
	start = init;
	end = endLe + 1;

	for (gen=start; gen<end; gen++)
	{
		if (gen<genLe)	expns[gen] = 1;
		else		expns[gen] = 0;

		if (gen==end-1)	getfiles[gen] = 1;
		else		getfiles[gen] = 0;
	}

	takeindividuals(gL1, motherL1, fatherL1, smatL1, qL1, initialgenL1);
}

/* ********************************************************************* */

MIXED_LINE()
{
	if (tracelevel!=0)	fprintf (fptr,"\n******* EXPERIMENT C (MIXED LINE) *******\n");

	bott = bottL;
	NIND = ib * bott;
	start = init;
	end = endmL + 1;

	for (gen=start; gen<end; gen++)
	{
		expns[gen] = 0;
		if (gen==end-1)	getfiles[gen] = 1;
		else		getfiles[gen] = 0;
	}

	mixLines();
}

/* ********************************************************************* */

sample_NP ()
{
//COMMENT
	if (tracelevel!=0)    fprintf(fptr,"\nCOMMENT - sampleNP\n");

	int np;
	
	/* ***** sample individuals from the Natural Population ***** */

	for (i=0; i<NIND; i++)
	{
		if(i%2==0)	do {ran_i = (int)(uniform() * NINDNP);}		while (ran_i%2!=0); /* sample females */
		else		do {ran_i = (int)(uniform() * NINDNP);}		while (ran_i%2==0); /* sample males */
		if (tracelevel!=0)	fprintf(fptr,"i=%d  ran_i=%d\n", i, ran_i);

		for (k=0; k<NCRO; k++)
		{
			np=gNP[i][k][0]; gNP[i][k][0]=gNP[ran_i][k][0]; gNP[ran_i][k][0]=np;
			np=gNP[i][k][1]; gNP[i][k][1]=gNP[ran_i][k][1]; gNP[ran_i][k][1]=np;
		}
	}

	for (i=0; i<NIND; i++)
	{
		for (k=0; k<NCRO; k++)
		{
			g[i][k][0]=gNP[i][k][0];
			g[i][k][1]=gNP[i][k][1];
		}
//COMMENT
//		if (tracelevel!=0)	fprintf(fptr,"g[%d][0][0]=%d\tg[%d][0][1]=%d\n", i, g[i][0][0], i, g[i][0][1]);
	}

	/* ********************* frequences ********************* */

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	{
		AA=0.0; Aa=0.0; aa=0.0;

		for (i=0; i<NIND; i++)
		{
			if (((g[i][k][0] & RM[l])==RM[l])&&((g[i][k][1] & RM[l])==RM[l]))	aa+=1.0;
			else if (((g[i][k][0] & RM[l])!=RM[l])&&((g[i][k][1] & RM[l])!=RM[l]))	AA+=1.0;
		     	else	Aa+=1.0;
		}

		q[k][l] = (aa/(double)NIND)+(Aa/(2.0*(double)NIND));

		initialgen[k][l] = initialgenNP[k][l];

//COMMENT
//		if (tracelevel!=0)	fprintf(fptr,"k=%d\tl=%d\tAA=%1.0f\tAa=%1.0f\taa=%1.0f\tq[k][l]=%f\tinitialgen[k][l]=%d\n",k,l,AA,Aa,aa,q[k][l],initialgen[k][l]);

	}
}

/* ********************************************************************* */

save(int sgnt[][MM][2], int *smoth, int *sfath, double ssmt[][NN], double qs[][31], int sinitgen[][31])
{
//COMMENT
//	if (tracelevel!=0)	fprintf(fptr,"\n\nCOMMENT - saveind\n");

	for (i=0; i<NIND; i++)
	{
		for (k=0; k<NCRO; k++)
		{
			sgnt[i][k][0]=g[i][k][0];
			sgnt[i][k][1]=g[i][k][1];
		}
		smoth[i]=mother[i];
		sfath[i]=father[i];
	}

	for (i=0; i<NIND; i++)
	for (j=0; j<NIND; j++)
	{
		ssmt[sfath[i]][sfath[j]]=smat[father[i]][father[j]];
		ssmt[sfath[i]][smoth[j]]=smat[father[i]][mother[j]];
		ssmt[smoth[i]][sfath[j]]=smat[mother[i]][father[j]];
		ssmt[smoth[i]][smoth[j]]=smat[mother[i]][mother[j]];
	}

	for (k=0; k<NCRO; k++) 
	for (l=0; l<NLOCI; l++)
	{
		qs[k][l] = q[k][l];
		sinitgen[k][l] = initialgen[k][l];
	}

//COMMENT
//	if (tracelevel!=0)	for (bt=0; bt<bott; bt++)	for (i=bt*ib; i<ib*(bt+1); i++)	fprintf(fptr,"bottle%d i=%d mother=%d father=%d  sgnt[%d][0][0]=%d  sgnt[%d][0][1]=%d  F=%f\n", bt, i, mother[i], father[i],  i, sgnt[i][0][0], i, sgnt[i][0][1], 2.0*(0.5 * (1.0 + ssmt[sfath[i]][smoth[i]]))-1.0); 
//	if (tracelevel!=0)	for (k=0; k<NCRO; k++)	for (l=0; l<NLOCI; l++)	fprintf(fptr,"k=%d l=%d qs[%d][%d]=%f sinitgen=%d\n", k, l, k, l, qs[k][l], sinitgen[k][l]);
}

/* ********************************************************************* */

sample(int btt, int tgnt[][MM][2], int *tmoth, int *tfath, double tsmt[][NN], double tq[][31], int tinitgen[][31])
{
	/* ******* SAMPLE RANDOM INDIVIDUALS, REGARDLESS THE BOTTLE ******* */

//COMMENT
	if (tracelevel!=0)	fprintf(fptr,"\n\nCOMMENT - sample\n");

	int b, sN;
	sN = ib * btt;

	for (i=0; i<NIND; i++)
	{
		if(i%2==0)	do {ran_i = (int)(uniform() * sN);}	while (ran_i%2!=0); /* sample females */
		else		do {ran_i = (int)(uniform() * sN);}	while (ran_i%2==0); /* sample males */
		if (tracelevel!=0)	fprintf(fptr,"i=%d  ran_i=%d\n", i, ran_i);

		for (k=0; k<NCRO; k++)
		{
			b=tgnt[i][k][0]; tgnt[i][k][0]=tgnt[ran_i][k][0]; tgnt[ran_i][k][0]=b;
			b=tgnt[i][k][1]; tgnt[i][k][1]=tgnt[ran_i][k][1]; tgnt[ran_i][k][1]=b;
		}

		b=tfath[i]; tfath[i]=tfath[ran_i]; tfath[ran_i]=b;
	   	b=tmoth[i]; tmoth[i]=tmoth[ran_i]; tmoth[ran_i]=b;
	}

	for (i=0; i<NIND; i++)
	{
		for (k=0; k<NCRO; k++)
		{
			g[i][k][0]=tgnt[i][k][0];
			g[i][k][1]=tgnt[i][k][1];
		}

		mother[i]=tmoth[i];
		father[i]=tfath[i];
	}

	for (i=0; i<NIND; i++)
	for (j=0; j<NIND; j++)
	{
		smat[father[i]][father[j]]=tsmt[tfath[i]][tfath[j]];
		smat[father[i]][mother[j]]=tsmt[tfath[i]][tmoth[j]];
		smat[mother[i]][father[j]]=tsmt[tmoth[i]][tfath[j]];
		smat[mother[i]][mother[j]]=tsmt[tmoth[i]][tmoth[j]];
	}

	for (k=0; k<NCRO; k++) 
	for (l=0; l<NLOCI; l++)
	{
		q[k][l] = tq[k][l];
		initialgen[k][l] = tinitgen[k][l];
	}

//COMMENT
//	if (tracelevel!=0)	for (bt=0; bt<bott; bt++)	for (i=bt*ib; i<ib*(bt+1); i++)	fprintf(fptr,"bottle%d i=%d mother=%d father=%d  g[%d][0][0]=%d  g[%d][0][1]=%d  F=%f\n", bt, i, mother[i], father[i], i, g[i][0][0], i, g[i][0][1], 2.0*(0.5 * (1.0 + smat[father[i]][mother[i]]))-1.0); 
//	if (tracelevel!=0)	for (k=0; k<NCRO; k++)		for (l=0; l<NLOCI; l++)	fprintf(fptr,"k=%d l=%l q[%d][%d]=%f initialgen=%d\n", k, l, k, l, q[k][l], initialgen[k][l]);

}

/* ********************************************************************* */

takeindividuals(int tgnt[][MM][2], int *tmoth, int *tfath, double tsmt[][NN], double tq[][31], int tinitgen[][31])
{
//COMMENT
//	if (tracelevel!=0)	fprintf(fptr,"\n\nCOMMENT - take saved individuals\n");

	for (i=0; i<NIND; i++)
	{
		for (k=0; k<NCRO; k++)
		{
			g[i][k][0]=tgnt[i][k][0];
			g[i][k][1]=tgnt[i][k][1];
		}

		mother[i]=tmoth[i];
		father[i]=tfath[i];
	}

	for (i=0; i<NIND; i++)
	for (j=0; j<NIND; j++)
	{
		smat[father[i]][father[j]]=tsmt[tfath[i]][tfath[j]];
		smat[father[i]][mother[j]]=tsmt[tfath[i]][tmoth[j]];
		smat[mother[i]][father[j]]=tsmt[tmoth[i]][tfath[j]];
		smat[mother[i]][mother[j]]=tsmt[tmoth[i]][tmoth[j]];
	}

	for (k=0; k<NCRO; k++) 
	for (l=0; l<NLOCI; l++)
	{
		q[k][l] = tq[k][l];
		initialgen[k][l] = tinitgen[k][l];
	}

//COMMENT
//	if (tracelevel!=0)	for (bt=0; bt<bott; bt++)	for (i=bt*ib; i<ib*(bt+1); i++)	fprintf(fptr,"bottle%d i=%d mother=%d father=%d  g[%d][0][0]=%d  g[%d][0][1]=%d  F=%f\n", bt, i, mother[i], father[i], i, g[i][0][0], i, g[i][0][1], 2.0*(0.5 * (1.0 + smat[father[i]][mother[i]]))-1.0); 
//	if (tracelevel!=0)	for (k=0; k<NCRO; k++)		for (l=0; l<NLOCI; l++)	fprintf(fptr,"k=%d l=%l q[%d][%d]=%f initialgen=%d\n", k, l, k, l, q[k][l], initialgen[k][l]);

}

/* ********************************************************************* */

mixLines()
{
//COMMENT
	if (tracelevel!=0)	fprintf(fptr,"COMMENT - mixLines\n");

	int li;

	/* ******** take individuals ******** */

	for (i=0; i<NIND; i++)
	{
		for (k=0; k<NCRO; k++)	
		{
			g[i][k][0]=gL1[i][k][0];
			g[i][k][1]=gL1[i][k][1]; 
		}	
		mother[i]=motherL1[i];
		father[i]=fatherL1[i];
	}
	for (i=NIND; i<NIND*2; i++)
	{
		for (k=0; k<NCRO; k++)
		{
			g[i][k][0]=gL2[i-NIND][k][0];
			g[i][k][1]=gL2[i-NIND][k][1];
		}
		mother[i]=motherL2[i-NIND] + NIND;
		father[i]=fatherL2[i-NIND] + NIND;
	}

	for (i=0; i<NIND; i++)
	for (j=0; j<NIND; j++)
	{
		smat[father[i]][father[j]]=smatL1[fatherL1[i]][fatherL1[j]];
		smat[father[i]][mother[j]]=smatL1[fatherL1[i]][motherL1[j]];
		smat[mother[i]][father[j]]=smatL1[motherL1[i]][fatherL1[j]];
		smat[mother[i]][mother[j]]=smatL1[motherL1[i]][motherL1[j]];
	}

	for (i=NIND; i<NIND*2; i++)
	for (j=NIND; j<NIND*2; j++)
	{
		smat[father[i]][father[j]]=smatL2[fatherL2[i-NIND]][fatherL2[j-NIND]];
		smat[father[i]][mother[j]]=smatL2[fatherL2[i-NIND]][motherL2[j-NIND]];
		smat[mother[i]][father[j]]=smatL2[motherL2[i-NIND]][fatherL2[j-NIND]];
		smat[mother[i]][mother[j]]=smatL2[motherL2[i-NIND]][motherL2[j-NIND]];
	}

//COMMENT
//	if (tracelevel!=0)	fprintf(fptr,"\nindividuals before randomize\n");
//	if (tracelevel!=0)	for (i=0; i<NIND; i++)	fprintf(fptr,"L1 i%d mother=%d father=%d - g[%d][0][0]=%d  g[%d][0][1]=%d\n", i, mother[i], father[i], i, g[i][0][0], i, g[i][0][1]);
//	if (tracelevel!=0)	for (i=NIND; i<NIND*2; i++)	fprintf(fptr,"L2 i%d mother=%d father=%d - g[%d][0][0]=%d  g[%d][0][1]=%d\n", i, mother[i], father[i], i, g[i][0][0], i, g[i][0][1]);

	/* ****** randomize and sample ****** */

	for (i=0; i<NIND; i++)
	{
		if(i%2==0)	do {ran_i = (int)(uniform() * (NIND*2));}	while (ran_i%2!=0); /* sample females */
		else		do {ran_i = (int)(uniform() * (NIND*2));}	while (ran_i%2==0); /* sample males */
		if (tracelevel!=0)	fprintf(fptr,"i=%d  ran_i=%d\n", i, ran_i);


		for (k=0; k<NCRO; k++)
		{
			li=g[i][k][0]; g[i][k][0]=g[ran_i][k][0]; g[ran_i][k][0]=li;
			li=g[i][k][1]; g[i][k][1]=g[ran_i][k][1]; g[ran_i][k][1]=li;
		}

		li=father[i]; father[i]=father[ran_i]; father[ran_i]=li;
		li=mother[i]; mother[i]=mother[ran_i]; mother[ran_i]=li;
	}

	/* *********** frequences *********** */

	for (k=0; k<NCRO; k++) 
	for (l=0; l<NLOCI; l++)
	{							
		AA=0.0; Aa=0.0; aa=0.0;
		for (i=0; i<NIND; i++)
		{
			if (((g[i][k][0] & RM[l])==RM[l])&&((g[i][k][1] & RM[l])==RM[l]))	aa+=1.0;
			else if (((g[i][k][0] & RM[l])!=RM[l])&&((g[i][k][1] & RM[l])!=RM[l]))	AA+=1.0;
			else	Aa+=1.0;
		}
		q[k][l] = (aa/((double)NIND))+(Aa/(2.0*(double)NIND));	

		if ( ((initialgenL1[k][l]<initialgenL2[k][l])&&(initialgenL1[k][l]==(-99))) || ((initialgenL1[k][l]>initialgenL2[k][l])&&(initialgenL2[k][l]!=(-99))) )		initialgen[k][l] = initialgenL2[k][l];
		else		initialgen[k][l] = initialgenL1[k][l];	
	}
	
//COMMENT
//	if (tracelevel!=0)	fprintf(fptr,"\nindividuals after sampling\n");
//	if (tracelevel!=0)	for (bt=0; bt<bott; bt++)	for (i=bt*ib; i<ib*(bt+1); i++)	fprintf(fptr,"bottle%d i%d mother=%d father=%d - g[%d][0][0]=%d  g[%d][0][1]=%d\n", bt, i, mother[i], father[i], i, g[i][0][0], i, g[i][0][1]);
//	if (tracelevel!=0)	for (k=0; k<NCRO; k++) 	for (l=0; l<NLOCI; l++)	fprintf(fptr,"k=%d l=%d q[%d][%d]=%f initialgen=%d\n", k, l, k, l, q[k][l], initialgen[k][l]);

}

/* ********************************************************************* */

poisson_tables ()
{   
	/* LOCI WITH POISSON (2NL) NEW MUTATIONS */

	if ( (exp(-2.0*(double)NIND*Lambda) != 0.0)&&(2.0*(double)NIND*Lambda < normalthreshold) )
	generatepoissontable(2.0*(double)NIND*Lambda, &lastinmutantspoissontable, mutantspoissontable, maxmpt-1);

	/* NUMERO DE RECOMBINACIONES POISSON CON MEDIA L*/

	if ( (exp(-L) != 0.0) && (L < normalthreshold) )
	generatepoissontable(L, &lastinrecombinantpoissontable, recombinantpoissontable, maxmpt-1);
}

/* ********************************************************************* */

mutation_neutral()
{
	/* NEUTRAL GENES: (POISSON) 2N(Lambda) NEW MUTATIONS (ONE DIRECTION) */

	muts = mutationnumber();

	if (tracelevel!=0)    fprintf(fptr,"\n New neutral mutants = %d\n", muts);

	countNoSS = 0;

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	{
		if (q[k][l]==0.0)
		{
			countNoSS += 1;
			NoSS_k[countNoSS-1] = k;
			NoSS_l[countNoSS-1] = l;
		}
	}
	if (tracelevel!=0)    fprintf(fptr," countNoSS=%d\n", countNoSS);

	mutants_ocurred = 0;

	if (countNoSS != 0)
	{
		disorder_NoSS (NoSS_k,NoSS_l);

		for (m=0; m<countNoSS; m++)
		{
			if (mutants_ocurred==muts)    goto label;

			ran_i = (int)(uniform()*NIND);
   			ran_h = (int)(uniform()*2.0);

			g[ran_i][NoSS_k[m]][ran_h]=(g[ran_i][NoSS_k[m]][ran_h] | RM[NoSS_l[m]]); 

		   	mutants_ocurred += 1;
			initialgen[NoSS_k[m]][NoSS_l[m]] = gen;
//COMMENT
//			if (tracelevel!=0)    fprintf(fptr,"ran_i=%d  k=%d  l=%d  ran_h=%d  initialgen=%d\n", ran_i, NoSS_k[m],  NoSS_l[m], ran_h, initialgen[NoSS_k[m]][NoSS_l[m]]);

		}
	}

	for (m=mutants_ocurred; m<muts; m++)
	{
		ran_i = (int)(uniform()*NIND);
		ran_k = (int)(uniform()*NCRO);
		ran_l = (int)(uniform()*NLOCI);
		ran_h = (int)(uniform()*2.0);

		g[ran_i][ran_k][ran_h]=(g[ran_i][ran_k][ran_h] | RM[ran_l]); 
//COMMENT
//		if (tracelevel!=0)    fprintf(fptr,"(rec. mut.) ran_i=%d  ran_k=%d  ran_l=%d  ran_h=%d\n", ran_i, ran_k, ran_l, ran_h);
	}

	label: /* end of mutations */;
}

/* ********************************************************************* */

int mutationnumber()
{
	int r;
	if ((2.0*(double)NIND*Lambda < normalthreshold) && (exp(-2.0*(double)NIND*Lambda) != 0.0) )
	{
		r = poisson(lastinmutantspoissontable, mutantspoissontable);
	}
	else r = (int)( normal(2.0*(double)NIND*Lambda, sqrt(2.0*(double)NIND*Lambda)) );
	return(r);
}

/* ********************************************************************* */

void disorder_NoSS (NoSS_k,NoSS_l)
int NoSS_k[], NoSS_l[];
{
	int a, b, rnd;
	
	for (i=0; i<countNoSS-1; i++)
	{
	   rnd=(int)(uniform()*(countNoSS-i));
	   a=NoSS_k[countNoSS-1-i]; NoSS_k[countNoSS-1-i]=NoSS_k[rnd]; NoSS_k[rnd]=a;
	   b=NoSS_l[countNoSS-1-i]; NoSS_l[countNoSS-1-i]=NoSS_l[rnd]; NoSS_l[rnd]=b;
	}
}

/* ********************************************************************* */

neutral_genes()
{
	double AAbt[bott], Aabt[bott], aabt[bott], qbt[MM][31], meanQbt[MM][31], Hbt[bott], HObt[bott]; 
	double Q, Q2, Qbt[bott], Q2bt[bott],seg, HO, Fstjunk, HTjunk, Fstjunklines, AF, AFbottle;
	int genFst;

//COMMENT
	if (tracelevel!=0)    fprintf(fptr,"\nCOMMENT - Neutral genes\n");

	H = 0.0; HO = 0.0; HTjunk = 0.0; Q = 0.0; Q2 = 0.0; seg = 0.0;
	for (bt=0; bt<bott; bt++)	{ Hbt[bt] = 0.0; HObt[bt] = 0.0; Qbt[bt] = 0.0; Q2bt[bt] = 0.0; }
	Hw_btmean[lin][gen] = 0.0;
	HO_btmean[lin][gen] = 0.0;
	Vq_btmean[lin][gen] = 0.0;
	Qbottle[lin][gen] = 0.0;
	sumq[lin][gen] = 0.0;
	sumq2[lin][gen] = 0.0;

	if (lin == 1)	
	{
		HTjunklines[gen] = 0.0;
		for (k=0; k<NCRO; k++) for (l=0; l<NLOCI; l++)	meanQlines[gen][k][l] = 0.0; 
	}
	for (k=0; k<NCRO; k++) for (l=0; l<NLOCI; l++)	meanQbt[k][l] = 0.0;


	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	{
		if(initialgen[k][l] != (-99))
		{
			AA=0.0; Aa=0.0; aa=0.0;
			for (bt=0; bt<bott; bt++) { AAbt[bt]=0.0; Aabt[bt]=0.0; aabt[bt]=0.0; }

			for (bt=0; bt<bott; bt++)	
			for (i=bt*ib; i<ib*(bt+1); i++)
	    		{
				if (((g[i][k][0] & RM[l])==RM[l])&&((g[i][k][1] & RM[l])==RM[l]))	{ aa+=1.0; aabt[bt]+=1.0; }
				else if (((g[i][k][0] & RM[l])!=RM[l])&&((g[i][k][1] & RM[l])!=RM[l]))	{ AA+=1.0; AAbt[bt]+=1.0; }
	 			else	{ Aa+=1.0; Aabt[bt]+=1.0; }
	    		}

			/* global */
			q[k][l] = (aa/(double)NIND)+(Aa/(2.0*(double)NIND));
			H += 2.0 * q[k][l] * (1.0 - q[k][l]);
			HO += (Aa/(double)NIND);

			Q += q[k][l];
			Q2 += q[k][l]*q[k][l];
			if (q[k][l] > 0.0)	seg += 1.0;

			/* per bottle */
			for (bt=0; bt<bott; bt++)
			{
				qbt[k][l] = (aabt[bt]/(double)ib)+(Aabt[bt]/(2.0*(double)ib));
				Hbt[bt] += 2.0 * qbt[k][l] * (1.0 - qbt[k][l]);
				HObt[bt] += (Aabt[bt]/(double)ib);
				meanQbt[k][l] += qbt[k][l];  /* qbar in purgingv12.c */
				Qbt[bt] += qbt[k][l];
				Q2bt[bt] += qbt[k][l]*qbt[k][l];
			}

			/* between lines */
			if ((experiment == 99) && ((lin == 1) || (lin == 2)))	meanQlines[gen][k][l] += q[k][l];
			if ((tracelevel!=0) && ((experiment == 99) && ((lin == 1) || (lin == 2))))    fprintf(fptr,"gen=%d  L%d  q[%d][%d]=%f  H=%f\n", gen, lin, k, l, q[k][l], H);
	
//COMMENT
			if ((tracelevel!=0)&&(lin==0)&&(gen<=1))	fprintf(fptr,"k=%d\tl=%d\tAA=%1.0f\tAa=%1.0f\taa=%1.0f\tq[k][l]=%f\tinitialgen[k][l]=%d\n",k,l,AA,Aa,aa,q[k][l],initialgen[k][l]);
		}
	}

	/* ***** GLOBAL ***** */

	accum(&HOglobal[lin][gen], HO/(double)TOTLOCI);
	Hw[lin][gen] = H/(double)TOTLOCI;
	accum(&Hwglobal[lin][gen], H/(double)TOTLOCI);
//	accum(&NeH[lin][gen], 1.0 / ( 2.0 * (1.0 - pow((Hw[lin][gen]) / (Hw[0][0]), 1.0/gen))));
	accum(&NeH[lin][gen], 1.0 / ( 2.0 * (1.0 - ((Hw[lin][gen]) / (Hw[lin][gen-1]))) ));

	Qmean[lin][gen] = Q/(double)TOTLOCI;
	accum(&Qglobal[lin][gen], Q/(double)TOTLOCI);
	accum(&lociseg[lin][gen], seg);

	if (Q2==(Q*Q)/(double)TOTLOCI)	Vq[lin][gen] = 0.0;
	else	Vq[lin][gen] = (Q2-((Q*Q)/(double)TOTLOCI))/((double)TOTLOCI-1);
	accum(&Varq[lin][gen], Vq[lin][gen]);
//	accum(&Nevq[lin][gen], 1.0 / ( 2.0 * (1.0 - pow((1-(4*Vq[lin][gen])) / (1-(4*Vq[0][0])), 1.0/gen))));
//	AF = (Vq[lin][gen] - Vq[0][0]) / (gen*((Qmean[0][0]*(1-Qmean[0][0])) -  Vq[0][0]));
	AF = (Vq[lin][gen] - Vq[lin][gen-1]) / ((Qmean[0][0]*(1-Qmean[0][0])) -  Vq[lin][gen-1]);
	accum(&Nevq[lin][gen], 1/(2*AF));
	

	/* *** PER BOTTLE *** */

	for (bt=0; bt<bott; bt++)	Hw_btmean[lin][gen] += (Hbt[bt]/(double)TOTLOCI)/(double)bott;
	accum(&Hwbottle[lin][gen], Hw_btmean[lin][gen]);
//	accum(&NeH_bottle[lin][gen], 1.0 / ( 2.0 * (1.0 - pow((Hw_btmean[lin][gen]) / (Hw_btmean[0][0]), 1.0/gen))));
	accum(&NeH_bottle[lin][gen], 1.0 / ( 2.0 * (1.0 - ((Hw_btmean[lin][gen]) / (Hw_btmean[lin][gen-1]))) ));

	for (bt=0; bt<bott; bt++)	HO_btmean[lin][gen] += (HObt[bt]/(double)TOTLOCI)/(double)bott;
	accum(&HObottle[lin][gen], HO_btmean[lin][gen]);

	for (bt=0; bt<bott; bt++)	Qbottle[lin][gen] += (Qbt[bt]/(double)TOTLOCI)/(double)bott;
	accum(&meanQbottle[lin][gen], Qbottle[lin][gen]);

 	for (bt=0; bt<bott; bt++)
	{
		if (Q2bt[bt]==(Qbt[bt]*Qbt[bt])/(double)TOTLOCI)	Vq_btmean[lin][gen] += 0.0;
		else	Vq_btmean[lin][gen] += ((Q2bt[bt]-((Qbt[bt]*Qbt[bt])/(double)TOTLOCI))/((double)TOTLOCI-1)) / (double)bott;
	}
	accum(&Varq_bt[lin][gen], Vq_btmean[lin][gen]);
//	accum(&Nevq_bottle[lin][gen], 1.0 / ( 2.0 * (1.0 - pow((1-(4*Vq_btmean[lin][gen])) / (1-(4*Vq_btmean[0][0])), 1.0/gen))));
//	AFbottle = (Vq_btmean[lin][gen] - Vq_btmean[0][0]) / (gen*((Qbottle[0][0]*(1-Qbottle[0][0])) -  Vq_btmean[0][0]));
	AFbottle = (Vq_btmean[lin][gen] - Vq_btmean[lin][gen-1]) / ((Qbottle[0][0]*(1-Qbottle[0][0])) -  Vq_btmean[lin][gen-1]);
	accum(&Nevq_bottle[lin][gen], 1/(2*AFbottle));


	/* ******* FST ******* */
	
	if ((experiment == 0.1) && (gen >= tLINES) && (isolation == 0))	genFst = gen - tLINES; /* Bp bottles maintained isolated since tLINES */
	else	genFst = gen;

	for (k=0; k<NCRO; k++) for (l=0; l<NLOCI; l++)	HTjunk += 2.0 * (meanQbt[k][l]/bott) * ( 1.0 - (meanQbt[k][l]/bott));
	HTjk[lin][gen] = HTjunk/TOTLOCI;
	accum(&HT[lin][gen], HTjunk/TOTLOCI);
	Fstjunk = (HTjk[lin][gen] - Hw_btmean[lin][gen]) / HTjk[lin][gen];
	accum(&Fst[lin][gen], Fstjunk);
	accum(&NeFst[lin][gen], 1.0 / ( 2.0 * (1.0 - pow((1-Fstjunk), 1.0/genFst))));

	/* ******* FST between lines (if experiment = 99) ******* */

	if (lin == 2)
	{
		for (k=0; k<NCRO; k++) for (l=0; l<NLOCI; l++)	HTjunklines[gen] += 2.0 * (meanQlines[gen][k][l]/(double)lines) * ( 1.0 - (meanQlines[gen][k][l]/(double)lines));
		accum(&HTlines[gen], HTjunklines[gen]/TOTLOCI);
		Fstjunklines = ((HTjunklines[gen]/TOTLOCI) - ((Hw[1][gen]+Hw[2][gen])/(double)lines)) / (HTjunklines[gen]/TOTLOCI);
		accum(&Fst_lines[gen], Fstjunklines);
		accum(&NeFst_lines[gen], 1.0 / ( 2.0 * (1.0 - pow((1-Fstjunklines), 1.0/(gen-tLINES)))));

		if (tracelevel!=0)    fprintf(fptr,"gen=%d  L%d  HTjunklines=%f\n", gen, lin, HTjunklines[gen]);
	}

//COMMENT
	if (tracelevel!=0)    fprintf(fptr,"HO=%.10f  HObottle=%.10f  Hw=%.10f  NeH=%.10f  Hwbottle=%.10f  HT=%.10f  NeH_bottle=%.10f  Fst=%.10f  NeFst=%.10f  vq=%f  Nevq=%f  vq_bottle=%f  Nevq_bottle=%f  q=%f  qbt=%f\n", HO/(double)TOTLOCI, HO_btmean[lin][gen], H/(double)TOTLOCI, 1.0 / ( 2.0 * (1.0 - pow((Hw[lin][gen]) / (Hw[0][0]), 1.0/gen))), Hw_btmean[lin][gen], HTjunk/TOTLOCI, 1.0 / ( 2.0 * (1.0 - pow((Hw_btmean[lin][gen]) / (Hw_btmean[0][0]), 1.0/gen))), Fstjunk, 1.0 / ( 2.0 * (1.0 - pow((1-Fstjunk), 1.0/gen))), Vq[lin][gen], 1.0 / ( 2.0 * (1.0 - pow((1-(4*Vq[lin][gen])) / (1-(4*Vq[0][0])), 1.0/gen))), Vq_btmean[lin][gen], 1.0 / ( 2.0 * (1.0 - pow((1-(4*Vq_btmean[lin][gen])) / (1-(4*Vq_btmean[0][0])), 1.0/gen))), Qmean[lin][gen], Qbottle[lin][gen]);
	if (tracelevel!=0)    for (bt=0; bt<bott; bt++)	fprintf(fptr,"HO[bt%d]=%.10f  Hw[bt%d]=%.10f  vq[bt%d]=%f  qbt[bt%d]=%f\n", bt, HObt[bt]/(double)TOTLOCI, bt, Hbt[bt]/(double)TOTLOCI, bt, ((Q2bt[bt]-((Qbt[bt]*Qbt[bt])/(double)TOTLOCI))/((double)TOTLOCI-1)), bt, Qbt[bt]/(double)TOTLOCI);
}

/* ********************************************************************* */

dumpoffspringaftermutation()
{
	if (tracelevel==0)   return (0);

	fprintf(fptr,"\n Offspring after mutation \n");	

	for (bt=0; bt<bott; bt++)	for (i=bt*ib; i<ib*(bt+1); i++)	fprintf(fptr,"bottle%d (g0 g1)	%d	%d\n", bt, g[i][0][0], g[i][0][1]);
}

/* ********************************************************************* */

distribution_q()
{
//COMMENT
//	if (tracelevel!=0)	fprintf(fptr,"\nCOMMENT - distribution of q values\n");

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	if (initialgen[k][l] != (-99)) 
	{
		if (q[k][l] == 0.0)				accum (&q_00[lin][gen/tgen],  1.0);
		else if ( (q[k][l] > 0.0) && (q[k][l] < 0.1) )	accum (&q_00_01[lin][gen/tgen],  1.0);
		else if ( (q[k][l] >= 0.1) && (q[k][l] < 0.2) )	accum (&q_01_02[lin][gen/tgen],  1.0);
		else if ( (q[k][l] >= 0.2) && (q[k][l] < 0.3) )	accum (&q_02_03[lin][gen/tgen],  1.0);
		else if ( (q[k][l] >= 0.3) && (q[k][l] < 0.4) )	accum (&q_03_04[lin][gen/tgen],  1.0);
		else if ( (q[k][l] >= 0.4) && (q[k][l] < 0.5) )	accum (&q_04_05[lin][gen/tgen],  1.0);
		else if ( (q[k][l] >= 0.5) && (q[k][l] < 0.6) )	accum (&q_05_06[lin][gen/tgen],  1.0);
		else if ( (q[k][l] >= 0.6) && (q[k][l] < 0.7) )	accum (&q_06_07[lin][gen/tgen],  1.0);
		else if ( (q[k][l] >= 0.7) && (q[k][l] < 0.8) )	accum (&q_07_08[lin][gen/tgen],  1.0);
		else if ( (q[k][l] >= 0.8) && (q[k][l] < 0.9) )	accum (&q_08_09[lin][gen/tgen],  1.0);
		else if ( (q[k][l] >= 0.9) && (q[k][l] < 1.0) )	accum (&q_09_10[lin][gen/tgen],  1.0);
		else if (q[k][l] == 1.0)			accum (&q_10[lin][gen/tgen],  1.0);

//COMMENT
//		if (tracelevel!=0)	fprintf(fptr,"k=%d l=%d q=%f ini=%d\n", k, l, q[k][l], initialgen[k][l]);
	}
}

/* ********************************************************************* */

coancestry_matrix()
{
	if (gen==0)
	{
		for (i=0; i<NIND; i++)
		for (j=0; j<NIND; j++)
		{
			if (i == j)	mat[i][j] = 0.5;
			else		mat[i][j] = 0.0;
		}
	}
	else
	{
		for (i=0; i<NIND; i++)
		for (j=0; j<NIND; j++)
		{
			if (i == j)	mat[i][j] = 0.5 * (1.0 + smat[father[i]][mother[i]]);
			else		mat[i][j] = 0.25 * (smat[father[i]][father[j]] + smat[father[i]][mother[j]] + smat[mother[i]][father[j]] + smat[mother[i]][mother[j]]);
		}
	}

	for (i=0; i<NIND; i++)
	for (j=0; j<NIND; j++)
	{
		smat[i][j] = mat[i][j];

		accum(&fp[lin][gen], mat[i][j]); 
		accum(&Fp[lin][gen], (2.0*mat[i][i]-1.0)); 
	}

//COMMENT
//	if (tracelevel!=0)
//	{
//		fprintf(fptr, "\nCOMMENT - coancestry_matrix\n");   
//		for (bt=0; bt<bott; bt++)
//		for (i=bt*ib; i<ib*(bt+1); i++)
//		{
//			fprintf(fptr, "gen=%d  bottle%d  i=%d   fat=%d   mot=%d   Fp=%f\n",gen,bt,i,father[i],mother[i], 2.0*mat[i][i]-1.0);
//		}
//	}

}

/* ********************************************************************* */

dumpparents()
{
	if (tracelevel==0)   return (0);

	fprintf(fptr,"\nParents\n");
	for (bt=0; bt<bott; bt++)
	for (i=bt*ib; i<ib*(bt+1); i++)
	{
		fprintf(fptr,"\nbottle%d i=%d  mother=%d  father=%d  Fp=%f", bt, i, mother[i], father[i], 2.0*mat[i][i]-1.0);
//		fprintf(fptr,"  g[i][0][0]=%d g[i][0][1]=%d", g[i][0][0], g[i][0][1]);
	}
	fprintf(fptr,"\n");
}

/* ********************************************************************* */

data_files(FILE *fmap, FILE *fdata)
{
	int mark[NN], chr, pt, nsample, gmt, arm, armbreak;
	double q_df[MM][31], ran_g;

	// sample
	if ((CROX==0)&&(NIND<NINDs))	nsample = NIND;			/* autosomes - both females and males are sampled */
	else if ((CROX==1)&&((NIND/2)<NINDs))	nsample = NIND/2;	/* X chromosome - only females are sampled */
	else	nsample = NINDs;

	if (tracelevel!=0)	fprintf (fptr,"\n\ndata_files (sample n=%d)\n", nsample);

	for (i=0; i<NIND; i++)	mark[i]=0;
	for (i=0; i<nsample; i++)
	{
		if(CROX==0)	do {ran_i = (int)(uniform() * NIND);}	while (mark[ran_i]!=0);
		else		do {ran_i = (int)(uniform() * NIND);}	while ((mark[ran_i]!=0)||(ran_i%2!=0));

		for (k=0; k<NCRO; k++)
		{
			gdf[i][k][0]=g[ran_i][k][0];
			gdf[i][k][1]=g[ran_i][k][1];
		}

		mark[ran_i]=1;

		if (tracelevel!=0)	fprintf (fptr,"s=%d  ran_i=%d\n", i, ran_i);
	}

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	{
		AA=0.0; Aa=0.0; aa=0.0;

		for (i=0; i<nsample; i++)
	    	{
			if (((gdf[i][k][0] & RM[l])==RM[l])&&((gdf[i][k][1] & RM[l])==RM[l]))	aa+=1.0;
			else if (((gdf[i][k][0] & RM[l])!=RM[l])&&((gdf[i][k][1] & RM[l])!=RM[l]))	AA+=1.0;
	 		else	Aa+=1.0;
	    	}

		q_df[k][l] = (aa/(double)nsample)+(Aa/(2.0*(double)nsample));
		if (q_df[k][l] > 0.0)	SNP[k][l] = 1; 
		else 			SNP[k][l] = 0;
	}

	// data.map

	pt = (int) ((L * 100000000) / (KCRO*30)); if (tracelevel!=0) fprintf (fptr,"\n\npt=%d\n",pt);
	ss = 0;

	arm = 0;
	armbreak = KCRO/2.0;
	for (chr=0; chr<CRO; chr++)
	for (k=chr*KCRO; k<KCRO*(chr+1); k++)
	{
		/* separate chromosomes in arms */
		if ( (k == chr*KCRO) || (k == ((chr*KCRO) + armbreak)) )	arm += 1;

		for (l=0; l<NLOCI; l++)
		{
			ss ++;
			if ((CROX == 1) && (SNP[k][l] == 1))	fprintf(fmap,"%d\tSNP%d\t0\t%d\n", chr+1, (ss-1)*pt, ((k*30)+l)*pt);
			if ((CROX == 0) && (SNP[k][l] == 1))	fprintf(fmap,"%d\tSNP%d\t0\t%d\n", arm, (ss-1)*pt, ((k*30)+l)*pt);
		}
	}

	// data.ped

	for (i=0; i<nsample; i++)
	{
		fprintf(fdata,"%d IND%d 0 0 0 -9 ", i+1, i+1);
		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		if (SNP[k][l] == 1)
		{
			if (((gdf[i][k][0] & RM[l])==RM[l])&&((gdf[i][k][1] & RM[l])==RM[l]))		fprintf(fdata,"2 2 ");
	    		else if (((gdf[i][k][0] & RM[l])!=RM[l])&&((gdf[i][k][1] & RM[l])!=RM[l]))	fprintf(fdata,"1 1 ");
	    		else if (((gdf[i][k][0] & RM[l])==RM[l])&&((gdf[i][k][1] & RM[l])!=RM[l]))	fprintf(fdata,"2 1 ");
	    		else if (((gdf[i][k][0] & RM[l])!=RM[l])&&((gdf[i][k][1] & RM[l])==RM[l]))	fprintf(fdata,"1 2 ");
		}
		fprintf(fdata,"\n");
	}

	return(0);
}

/* ********************************************************************* */

mating()
{
	int mo, fa, EE[MM], FF[MM], family[NN];
	int numberrecs, nr, pointrec[MM][31], ncrorec[MM], rndk, rndl, marker, chr;
	int mark[NN], cm, ran_bt, fp, lp;
	double rnd;
	double family_sum=0.0, family_sum2=0.0;

	/* ************** KEEP GENOTYPES OF PARENTS ************** */ 

	for (i=0; i<NIND; i++)
	for (k=0; k<NCRO; k++)
	{
		sp[i][k][0]=g[i][k][0];
		sp[i][k][1]=g[i][k][1];
	}

	for (i=0; i<NIND; i++)   family[i]=0;

	/* ****************** ROUNDS OF MATING ***************** */

	if (tracelevel!=0)	fprintf (fptr,"\n\nRounds of random mating\n");

	NOFF = ib;
	if (expns[gen]==0)	bottOFF = bott;
	else
	{
		if (expns[gen+1]==0)	bottOFF = maxexp;
		else			bottOFF = bott + rate;
	}

	for (bt=0; bt<bottOFF; bt++)
	for (i=bt*NOFF; i<NOFF*(bt+1); i++)
	{
		/* ********** parents ********** */

		/* complete panmixia */
		if (matingsystem == 2) 
		{
			do {mo = (uniform() * NIND);} while (mo%2!=0); /* female */
			do {fa = (uniform() * NIND);} while (fa%2==0); /* male */

			mother[i] = mo;
			father[i] = fa;
			family[mo]+=1;
			family[fa]+=1;

			if (tracelevel!=0)	fprintf (fptr,"Progeny %d  bottle%d  mo = %d    fa = %d\n", i, bt, mo, fa);
		}
		else 
		{
			/* mating within bottles */
			if (bt<bott)
			{
				do {mo = (uniform() * NOFF) + (NOFF*bt);} while (mo%2!=0); /* female */
				do {fa = (uniform() * NOFF) + (NOFF*bt);} while (fa%2==0); /* male */

				mother[i] = mo;
				father[i] = fa;
				family[mo]+=1;
				family[fa]+=1;

				if (tracelevel!=0)	fprintf (fptr,"Progeny %d  bottle%d  mo = %d    fa = %d\n", i, bt, mo, fa);
			}
			else
			{
				ran_bt = (int)(uniform() * bott);

				do {mo = (int)(uniform() * NOFF) + (NOFF*ran_bt);} while (mo%2!=0); /* female */
				do {fa = (int)(uniform() * NOFF) + (NOFF*ran_bt);} while (fa%2==0); /* male */

				mother[i] = mo;
				father[i] = fa;
				family[mo]+=1;
				family[fa]+=1;

				if (tracelevel!=0)	fprintf (fptr,"Progeny %d  bottle%d mo = %d    fa = %d\n", i, bt, mo, fa);
			}
		}

		/* ******* Free recombination ******* */

		if(L==99.0)
		{
			if(CROX==0)
			{
				for (k=0; k<NCRO; k++)
				{
				   	EE[k] = (int)(uniform()*(pow(2.0,(double)NLOCI)));
					FF[k] = ~EE[k];
				   	g[i][k][0]=((EE[k]&sp[mo][k][0])|(FF[k]&sp[mo][k][1]));
				}
// COMMENT
//				if (tracelevel!=0)   fprintf (fptr,"  i=%d EE[0]=%d EE[1]=%d EE[2]=%d sf00=%d sf01=%d sf10=%d sf11=%d sf20=%d sf21=%d g00=%d g10=%d g20=%d \n", i, EE[0], EE[1], EE[2], sp[mo][0][0], sp[mo][0][1], sp[mo][1][0], sp[mo][1][1], sp[mo][2][0], sp[mo][2][1], g[i][0][0], g[i][1][0], g[i][2][0]);

				for (k=0; k<NCRO; k++)
				{
			 	  	EE[k] = (int)(uniform()*(pow(2.0,(double)NLOCI)));
			   		FF[k] = ~EE[k];
			  	 	g[i][k][1]=((EE[k]&sp[fa][k][0])|(FF[k]&sp[fa][k][1]));
				}
// COMMENT
//				if (tracelevel!=0)   fprintf (fptr,"i=%d EE[0]=%d EE[1]=%d EE[2]=%d sm00=%d sm01=%d sm10=%d sm11=%d sm20=%d sm21=%d g00=%d g10=%d g20=%d \n", i, EE[0], EE[1], EE[2], sp[fa][0][0], sp[fa][0][1], sp[fa][1][0], sp[fa][1][1], sp[fa][2][0], sp[fa][2][1], g[i][0][1], g[i][1][1], g[i][2][1]);
			}
			else
			{
				if(i%2==0)
				{
					/* FEMALE OFFSPRING */
					for (k=0; k<NCRO; k++)
					{
						/* mothers gamete */
					   	EE[k] = (int)(uniform()*(pow(2.0,(double)NLOCI)));
						FF[k] = ~EE[k];
		   				g[i][k][0]=((EE[k]&sp[mo][k][0])|(FF[k]&sp[mo][k][1]));
					}
//					if (tracelevel!=0)   fprintf (fptr,"  i=%d EE[0]=%d EE[1]=%d EE[2]=%d sf00=%d sf01=%d sf10=%d sf11=%d sf20=%d sf21=%d g00=%d g10=%d g20=%d \n", i, EE[0], EE[1], EE[2], sp[mo][0][0], sp[mo][0][1], sp[mo][1][0], sp[mo][1][1], sp[mo][2][0], sp[mo][2][1], g[i][0][0], g[i][1][0], g[i][2][0]);

					for (k=0; k<NCRO; k++)
					{
						/* fathers gamete */
		   				g[i][k][1]=sp[fa][k][0];
					}
				}
				else
				{
					/* MALE OFFSPRING */ 
					for (k=0; k<NCRO; k++)
					{
						/* mothers gamete */
					   	EE[k] = (int)(uniform()*(pow(2.0,(double)NLOCI)));
						FF[k] = ~EE[k];
		   				g[i][k][0]=((EE[k]&sp[mo][k][0])|(FF[k]&sp[mo][k][1]));

						/* no fathers gamete (X chromosome): g[i][k][1] = mothers gamete for simplicity */
						g[i][k][1]=((EE[k]&sp[mo][k][0])|(FF[k]&sp[mo][k][1]));
					}
//					if (tracelevel!=0)   fprintf (fptr,"  i=%d EE[0]=%d EE[1]=%d EE[2]=%d sf00=%d sf01=%d sf10=%d sf11=%d sf20=%d sf21=%d g00=%d g10=%d g20=%d \n", i, EE[0], EE[1], EE[2], sp[mo][0][0], sp[mo][0][1], sp[mo][1][0], sp[mo][1][1], sp[mo][2][0], sp[mo][2][1], g[i][0][0], g[i][1][0], g[i][2][0]);
				}
			}
		}

		/* **** Restricted recombination **** */ 

		else
		{
			/* Chromosome from mother */

			for (k=0; k<NCRO; k++)
			{
				ncrorec[k] = 0;
				for (l=0; l<NLOCI; l++)  pointrec[k][l] = 0;
			}

			for (chr=0; chr<CRO; chr++)
			{
				numberrecs = recombinationnumber();	
//				if (tracelevel!=0)   fprintf (fptr,"CRO%d\tL=%f  numberrecs (mo) = %d  ",chr,L,numberrecs);

				for (nr=0; nr<numberrecs; nr++)
				{
					rndk = (int)(uniform()*KCRO) + (KCRO*chr);
					rndl = (int)(uniform()*NLOCI);
					ncrorec[rndk] = 1;
					pointrec[rndk][rndl] = 1;
				}
			}
//			if (tracelevel!=0)   fprintf (fptr,"\n");

			marker = 1;

			for (k=0; k<NCRO; k++)
			{
				EE[k]=0;
				if (ncrorec[k] == 0)
				{
					if (marker==(-1))
					{
						EE[k] = ~EE[k];
					}
				}
				else
				{
					for (l=0; l<NLOCI; l++)
			      		{
						if (pointrec[k][l] == 0)
						{
							if (marker==(-1))  EE[k] = EE[k] | RM[l];
						}
						else
						{
							if (marker==1)
							{
								EE[k] = EE[k] | RM[l];
								marker = marker * (-1);
							}
							else
							{
								marker = marker * (-1);
							}
						}
					}
				}
				if ((k+1)%KCRO == 0)	if (uniform() < 0.5)	marker = marker * (-1);
			}

			for (chr=0; chr<CRO; chr++)
			{
				rnd = uniform();
				fp = chr*KCRO;
				lp = KCRO*(chr+1);

				for (k=fp; k<lp; k++)
				{
					/* mothers gamete */
					if (rnd < 0.5)	EE[k] = ~EE[k];
					FF[k] = ~EE[k];
					g[i][k][0]=((EE[k]&sp[mo][k][0])|(FF[k]&sp[mo][k][1]));

					if ((CROX==1)&&(i%2!=0))
					{
						/* male offspring so no fathers gamete (X chromosome): g[i][k][1] = mothers gamete for simplicity */
						g[i][k][1]=((EE[k]&sp[mo][k][0])|(FF[k]&sp[mo][k][1]));
					}
				}
			}
// COMMENT
/*			if (tracelevel!=0)   fprintf (fptr,"i=%d EE[0]=%d EE[1]=%d EE[2]=%d sf00=%d sf01=%d sf10=%d sf11=%d sf20=%d sf21=%d g00=%d g10=%d g20=%d \n", i, EE[0], EE[1], EE[2], sp[mo][0][0], sp[mo][0][1], sp[mo][1][0], sp[mo][1][1], sp[mo][2][0], sp[mo][2][1], g[i][0][0], g[i][1][0], g[i][2][0]);

			if (tracelevel!=0)
			{
				fprintf (fptr,"EE\n");
				for (k=0; k<NCRO; k++)
				{
					for (l=0; l<NLOCI; l++)
					{
						if((EE[k]&RM[l])==RM[l])  	fprintf (fptr,"1 ");
						else			   	fprintf (fptr,"0 ");
					}
					fprintf (fptr,"\n");
				}
				fprintf (fptr,"FF\n");
				for (k=0; k<NCRO; k++)
				{
					for (l=0; l<NLOCI; l++)
					{
						if((FF[k]&RM[l])==RM[l])  	fprintf (fptr,"1 ");
						else			   	fprintf (fptr,"0 ");
					}
					fprintf (fptr,"\n");
				}
			}
*/
			/* Chromosome from father */

			for (k=0; k<NCRO; k++)
			{
				ncrorec[k] = 0;
				for (l=0; l<NLOCI; l++)  pointrec[k][l] = 0;
			}

			for (chr=0; chr<CRO; chr++)
			{
				numberrecs = recombinationnumber();	
//				if (tracelevel!=0)   fprintf (fptr,"CRO%d\tL=%f  numberrecs (fa) = %d  ",chr,L,numberrecs);

				for (nr=0; nr<numberrecs; nr++)
				{
					rndk = (int)(uniform()*KCRO) + (KCRO*chr);
					rndl = (int)(uniform()*NLOCI);
					ncrorec[rndk] = 1;
					pointrec[rndk][rndl] = 1;
				}
			}
//			if (tracelevel!=0)   fprintf (fptr,"\n");

			marker = 1;

			for (k=0; k<NCRO; k++)
			{
				EE[k]=0;
				if (ncrorec[k] == 0)
				{
					if (marker==(-1))
					{
						EE[k] = ~EE[k];
					}
				}
				else
				{
					for (l=0; l<NLOCI; l++)
			      		{
						if (pointrec[k][l] == 0)
						{
							if (marker==(-1))  EE[k] = EE[k] | RM[l];
						}
						else
						{
							if (marker==1)
							{
								EE[k] = EE[k] | RM[l];
								marker = marker * (-1);
							}
							else
							{
								marker = marker * (-1);
							}
						}
					}
				}
				if ((k+1)%KCRO == 0)	if (uniform() < 0.5)	marker = marker * (-1);
			}

			for (chr=0; chr<CRO; chr++)
			{
				rnd = uniform();
				fp = chr*KCRO;
				lp = KCRO*(chr+1);

				for (k=fp; k<lp; k++)
				{
					/* fathers gamete */
					if (rnd < 0.5)	EE[k] = ~EE[k];
					FF[k] = ~EE[k];
					if (CROX==0)	g[i][k][1]=((EE[k]&sp[fa][k][0])|(FF[k]&sp[fa][k][1]));
					else
					{
						if(i%2==0)	g[i][k][1]=sp[fa][k][0];
						else	/* male offspring so no fathers gamete (X chromosome) */;
					}
				}
			}
// COMMENT
//			if (tracelevel!=0)	for (k=0; k<5; k++)	fprintf (fptr,"i%d k%d g[i][k][0]=%d  g[i][k][1]=%d   mo: sp[mo][k][0]=%d  sp[mo][k][1]=%d   fa: sp[fa][k][0]=%d  sp[fa][k][1]=%d\n", i, k, g[i][k][0], g[i][k][1], sp[mo][k][0], sp[mo][k][1], sp[fa][k][0], sp[fa][k][1]);
/*			if (tracelevel!=0)	fprintf (fptr,"i=%d EE[0]=%d EE[1]=%d EE[2]=%d sm00=%d sm01=%d sm10=%d sm11=%d sm20=%d sm21=%d g00=%d g10=%d g20=%d \n", i, EE[0], EE[1], EE[2], sp[fa][0][0], sp[fa][0][1], sp[fa][1][0], sp[fa][1][1], sp[fa][2][0], sp[fa][2][1], g[i][0][1], g[i][1][1], g[i][2][1]);
			if (tracelevel!=0)
			{
				fprintf (fptr,"EE\n");
				for (k=0; k<NCRO; k++)
				{
					for (l=0; l<NLOCI; l++)
					{
						if((EE[k]&RM[l])==RM[l])  	fprintf (fptr,"1 ");
						else			   	fprintf (fptr,"0 ");
					}
					fprintf (fptr,"\n");
				}
				fprintf (fptr,"FF\n");
				for (k=0; k<NCRO; k++)
				{
					for (l=0; l<NLOCI; l++)
					{
						if((FF[k]&RM[l])==RM[l])  	fprintf (fptr,"1 ");
						else			   	fprintf (fptr,"0 ");
					}
					fprintf (fptr,"\n");
				}
			}
*/
		}
	}

	if (matingsystem == 0)	circulardistribution();
	if (matingsystem == 1)	randomiseoffspring();

	/* ************** VARIANCE OF FAMILY SIZE ************** */

	for (i=0; i<NIND; i++)
	{
		family_sum += family[i];
		family_sum2 += (family[i] * family[i]);
	}

	sk2 = ((family_sum2-(family_sum*family_sum/NIND))/(NIND-1.0)); 

	accum(&SK2[lin][gen], sk2);

	ne = (4.0 * NIND) / (2.0 + sk2);
	accum(&Ne[lin][gen], ne);

	if (tracelevel!=0)	fprintf (fptr,"\nContributions from parents\n\n");
	if (tracelevel!=0)	for (i=0; i<NIND; i++)	fprintf (fptr,"family(%d) = %d\n", i, family[i]);
	if (tracelevel!=0)	fprintf (fptr,"\nSK2 = %5.3f\n", sk2);

	NIND = NOFF * bottOFF;
	bott = bottOFF;
}

/* ***************************************************** */

circulardistribution()
{
	/* ******** CIRCULAR DISTRIBUTION OF OFFSPRING ********* */

//COMMENT
	if (tracelevel!=0)	fprintf(fptr,"\nOffspring BEFORE circular distribution\n");
	if (tracelevel!=0)	for (i=0; i<(bottOFF*NOFF); i++)	fprintf(fptr,"i=%d  mother=%d  father=%d	g[i][0][0]=%d	g[i][0][1]=%d	F=%f\n", i, mother[i], father[i], g[i][0][0], g[i][0][1], 2.0*(0.5 * (1.0 + smat[father[i]][mother[i]]))-1.0);

// bott instead of bottOFF to exclude the duplicated bottles in the duplication generation, because they are already a pool of the others
	
	if ((experiment==0.1) && (gen >= tLINES))	mig = migration * isolation;
	else	mig = migration;

	for (bt=0; bt<bott; bt++)
	for (i=bt*NOFF; i<NOFF*(bt+1); i++)
	{
		if (i >= ((NOFF*bt)+(NOFF*mig)))
		{
			for (k=0; k<NCRO; k++)
			{
				sp[i][k][0]=g[i][k][0];
				sp[i][k][1]=g[i][k][1];
			}
			spmother[i] = mother[i];
			spfather[i] = father[i];
		}
		else
		{
			if (bt == (bott-1))
			{
				for (k=0; k<NCRO; k++)
				{
					sp[i][k][0]=g[i-(NOFF*bt)][k][0];
					sp[i][k][1]=g[i-(NOFF*bt)][k][1];
				}
				spmother[i] = mother[i-(NOFF*bt)];
				spfather[i] = father[i-(NOFF*bt)];
			}
			else
			{
				for (k=0; k<NCRO; k++)
				{
					sp[i][k][0]=g[i+NOFF][k][0];
					sp[i][k][1]=g[i+NOFF][k][1];
				}
				spmother[i] = mother[i+NOFF];
				spfather[i] = father[i+NOFF];
			}
		}
	}

	for (i=0; i<(bott*NOFF); i++)
	{
		for (k=0; k<NCRO; k++)
		{
			g[i][k][0]=sp[i][k][0];
			g[i][k][1]=sp[i][k][1];
		}
		mother[i] = spmother[i];
		father[i] = spfather[i];
	}

//COMMENT
	if (tracelevel!=0)	fprintf(fptr,"\nOffspring AFTER circular distribution\n");
	if (tracelevel!=0)	for (i=0; i<(bottOFF*NOFF); i++)	fprintf(fptr,"i=%d  mother=%d  father=%d	g[i][0][0]=%d	g[i][0][1]=%d	F=%f\n", i, mother[i], father[i], g[i][0][0], g[i][0][1], 2.0*(0.5 * (1.0 + smat[father[i]][mother[i]]))-1.0);
}

/* ********************************************************************* */

randomiseoffspring()
{
//COMMENT
	if (tracelevel!=0)	fprintf(fptr,"\nOffspring BEFORE randomization\n");
//	if (tracelevel!=0)	for (i=0; i<(bottOFF*NOFF); i++)	fprintf(fptr,"i=%d  mother=%d  father=%d	g[i][0][0]=%d	g[i][0][1]=%d	F=%f\n", i, mother[i], father[i], g[i][0][0], g[i][0][1], 2.0*(0.5 * (1.0 + smat[father[i]][mother[i]]))-1.0);

	/* ************ RANDOMISE PROGENY ************ */

	for (i=0; i<(bottOFF*NOFF); i++)
	{
		if(i%2==0)	do {ran_i=(int)(uniform()*(bottOFF*NOFF));}	while (ran_i%2!=0);	/* females */
		else		do {ran_i=(int)(uniform()*(bottOFF*NOFF));}	while (ran_i%2==0);	/* males */
		if (tracelevel!=0)	fprintf(fptr,"i=%d  ran_i=%d\n", i, ran_i);

		for (k=0; k<NCRO; k++)
		{
			a=g[i][k][0]; g[i][k][0]=g[ran_i][k][0]; g[ran_i][k][0]=a;
			a=g[i][k][1]; g[i][k][1]=g[ran_i][k][1]; g[ran_i][k][1]=a;
		}
		a=father[i]; father[i]=father[ran_i]; father[ran_i]=a;
		a=mother[i]; mother[i]=mother[ran_i]; mother[ran_i]=a;
	}

//COMMENT
//	if (tracelevel!=0)	fprintf(fptr,"\nOffspring AFTER randomization\n");
//	if (tracelevel!=0)	for (i=0; i<(bottOFF*NOFF); i++)	fprintf(fptr,"i=%d  mother=%d  father=%d	g[i][0][0]=%d	g[i][0][1]=%d	F=%f\n", i, mother[i], father[i], g[i][0][0], g[i][0][1], 2.0*(0.5 * (1.0 + smat[father[i]][mother[i]]))-1.0);
}

/* ********************************************************************* */

int recombinationnumber ()
{
	int r;
	if ((L < normalthreshold) && (exp(-L) != 0.0) )
	{
		r = poisson(lastinrecombinantpoissontable, recombinantpoissontable);
	}
	else r = (int)normal(L, sqrt(L));
	return(r);
}

/* ********************************************************************* */

dumpoffspring()
{
	if (tracelevel==0)   return (0);

	fprintf(fptr,"\n Offspring before mutation \n");

	for (i=0; i<NIND; i++)	fprintf(fptr,"(g0 g1)	%d	%d\n", g[i][0][0], g[i][0][1]);
}

/* ********************************************************************* */

printout()
{
	if ((experiment!=0.2)&&(experiment!=2)&&(experiment!=3)) {
		fprintf(fgenBP,"** Natural Population N=%d\n", NINDNP);
		fprintf(fgenBP,"** Base Population N=%d (%d bottles)\n",ib*bottBP,bottBP);
		fprintf(fgenBP,"** Mating system %d\n",matingsystem); }
	if ((experiment==0.2)||(experiment==2)||(experiment==99)) {
		fprintf(fgenL,"** Natural Population N=%d\n", NINDNP);
		fprintf(fgenL,"** Base Population N=%d (%d bottles)\n",ib*bottBP,bottBP); 
		fprintf(fgenL,"** Mating system %d\n",matingsystem); }		
	if ((experiment==3)||(experiment==99)) {		
		fprintf(fgenmL,"** Natural Population N=%d\n", NINDNP);
		fprintf(fgenmL,"** Base Population N=%d (%d bottles)\n",ib*bottBP,bottBP);
		fprintf(fgenmL,"** Formation of lines: tLINES = %d N = %d (%d bottles)\n",tLINES,ib*bottL,bottL); 
		fprintf(fgenmL,"** Mating system %d\n",matingsystem); }
		
	if ((experiment!=0.2)&&(experiment!=2)&&(experiment!=3)) 	fprintf(fgenBP,"L=%4.2f  N.N-LOCI=%d  NSEGLOCNP=%d  nCHROM=%d (cro.X=%d)  Lambda(neutral model)=%f  N_sample(data files)=%d  exchange=%f  reps=%d\n\n", L, TOTLOCI, NSEGLOCNP, CRO, CROX, Lambda, NINDs, mig, replicates);
	if ((experiment==0.2)||(experiment==2)||(experiment==99))	fprintf(fgenL,"L=%4.2f  N.N-LOCI=%d  NSEGLOCNP=%d  nCHROM=%d (cro.X=%d)  Lambda(neutral model)=%f  N_sample(data files)=%d  exchange=%f  reps=%d\n\n", L, TOTLOCI, NSEGLOCNP, CRO, CROX, Lambda, NINDs, mig, replicates);
	if ((experiment==3)||(experiment==99))				fprintf(fgenmL,"L=%4.2f  N.N-LOCI=%d  NSEGLOCNP=%d  nCHROM=%d (cro.X=%d)  Lambda(neutral model)=%f  N_sample(data files)=%d  exchange=%f  reps=%d\n\n", L, TOTLOCI, NSEGLOCNP, CRO, CROX, Lambda, NINDs, mig, replicates);	

	//EXPERIMENT A
	if (experiment==0.1)
	{
		fprintf(fgenBP,"Base Population\n");
		fprintf(fgenBP,"%-6s  %-15s  %-10s  %-15s  %-10s  %-10s  %-10s  %-15s  %-10s  %-10s  %-10s  %-10s  %-15s  %-10s  %-10s  %-10s  %-15s  %-10s  %-15s  %-15s  %-10s  %-15s\n", "gen", "N", "bottles", "NeVk", "fp", "Fp", "SK2", "loci_seg", "Qglobal", "Qbottle", "HO", "Hw", "Hw(bottle)", "HT", "Fst", "Varq", "Varq(bottle)", "NeH", "NeH(bottle)", "NeFst", "NeVq", "NeVq(bottle)");
		for (gen=0; gen<init; gen++)		fprintf(fgenBP,"%-6d  %-15f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-15f  %-15f  %-10f  %-15f\n", gen, accmean(&N[0][gen]), accmean(&nbott[0][gen]), accmean(&Ne[0][gen]), accmean(&fp[0][gen]), accmean(&Fp[0][gen]), accmean(&SK2[0][gen]), accmean(&lociseg[0][gen]), accmean(&Qglobal[0][gen]), accmean(&meanQbottle[0][gen]), accmean(&HOglobal[0][gen]), accmean(&Hwglobal[0][gen]), accmean(&Hwbottle[0][gen]), accmean(&HT[0][gen]), accmean(&Fst[0][gen]), accmean(&Varq[0][gen]), accmean(&Varq_bt[0][gen]), accmean(&NeH[0][gen]), accmean(&NeH_bottle[0][gen]), accmean(&NeFst[0][gen]), accmean(&Nevq[0][gen]), accmean(&Nevq_bottle[0][gen]));
		fprintf(fgenBP,"\n");
	}
	if (experiment==1.0)
	{
		fprintf(fgenBP,"Reduction of Base Population\n");
		fprintf(fgenBP,"%-6s  %-15s  %-10s  %-15s  %-10s  %-10s  %-10s  %-15s  %-10s  %-10s  %-10s  %-10s  %-15s  %-10s  %-10s  %-10s  %-15s  %-10s  %-15s  %-15s  %-10s  %-15s\n", "gen", "N", "bottles", "NeVk", "fp", "Fp", "SK2", "loci_seg", "Qglobal", "Qbottle", "HO", "Hw", "Hw(bottle)", "HT", "Fst", "Varq", "Varq(bottle)", "NeH", "NeH(bottle)", "NeFst", "NeVq", "NeVq(bottle)");
		for (gen=0; gen<init; gen++)		fprintf(fgenBP,"%-6d  %-15f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-15f  %-15f  %-10f  %-15f\n", gen, accmean(&N[0][gen]), accmean(&nbott[0][gen]), accmean(&Ne[0][gen]), accmean(&fp[0][gen]), accmean(&Fp[0][gen]), accmean(&SK2[0][gen]), accmean(&lociseg[0][gen]), accmean(&Qglobal[0][gen]), accmean(&meanQbottle[0][gen]), accmean(&HOglobal[0][gen]), accmean(&Hwglobal[0][gen]), accmean(&Hwbottle[0][gen]), accmean(&HT[0][gen]), accmean(&Fst[0][gen]), accmean(&Varq[0][gen]), accmean(&Varq_bt[0][gen]), accmean(&NeH[0][gen]), accmean(&NeH_bottle[0][gen]), accmean(&NeFst[0][gen]), accmean(&Nevq[0][gen]), accmean(&Nevq_bottle[0][gen]));
		for (gen=init; gen<genBPr; gen++)	fprintf(fgenBP,"%-6d  %-15f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-15f  %-15f  %-10f  %-15f\n", gen, accmean(&N[3][gen]), accmean(&nbott[3][gen]), accmean(&Ne[3][gen]), accmean(&fp[3][gen]), accmean(&Fp[3][gen]), accmean(&SK2[3][gen]), accmean(&lociseg[3][gen]), accmean(&Qglobal[3][gen]), accmean(&meanQbottle[3][gen]), accmean(&HOglobal[3][gen]), accmean(&Hwglobal[3][gen]), accmean(&Hwbottle[3][gen]), accmean(&HT[3][gen]), accmean(&Fst[3][gen]), accmean(&Varq[3][gen]), accmean(&Varq_bt[3][gen]), accmean(&NeH[3][gen]), accmean(&NeH_bottle[3][gen]), accmean(&NeFst[3][gen]), accmean(&Nevq[3][gen]), accmean(&Nevq_bottle[3][gen]));
		fprintf(fgenBP,"\n");
	}
	if ((experiment==1.1)||(experiment==99))
	{
		fprintf(fgenBP,"Bottleneck & Expansion (expA)\n");
		fprintf(fgenBP,"%-6s  %-15s  %-10s  %-15s  %-10s  %-10s  %-10s  %-15s  %-10s  %-10s  %-10s  %-10s  %-15s  %-10s  %-10s  %-10s  %-15s  %-10s  %-15s  %-15s  %-10s  %-15s\n", "gen", "N", "bottles", "NeVk", "fp", "Fp", "SK2", "loci_seg", "Qglobal", "Qbottle", "HO", "Hw", "Hw(bottle)", "HT", "Fst", "Varq", "Varq(bottle)", "NeH", "NeH(bottle)", "NeFst", "NeVq", "NeVq(bottle)");
		for (gen=0; gen<init; gen++)		fprintf(fgenBP,"%-6d  %-15f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-15f  %-15f  %-10f  %-15f\n", gen, accmean(&N[0][gen]), accmean(&nbott[0][gen]), accmean(&Ne[0][gen]), accmean(&fp[0][gen]), accmean(&Fp[0][gen]), accmean(&SK2[0][gen]), accmean(&lociseg[0][gen]), accmean(&Qglobal[0][gen]), accmean(&meanQbottle[0][gen]), accmean(&HOglobal[0][gen]), accmean(&Hwglobal[0][gen]), accmean(&Hwbottle[0][gen]), accmean(&HT[0][gen]), accmean(&Fst[0][gen]), accmean(&Varq[0][gen]), accmean(&Varq_bt[0][gen]), accmean(&NeH[0][gen]), accmean(&NeH_bottle[0][gen]), accmean(&NeFst[0][gen]), accmean(&Nevq[0][gen]), accmean(&Nevq_bottle[0][gen]));
		for (gen=init; gen<genBPb; gen++)	fprintf(fgenBP,"%-6d  %-15f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-15f  %-15f  %-10f  %-15f\n", gen, accmean(&N[3][gen]), accmean(&nbott[3][gen]), accmean(&Ne[3][gen]), accmean(&fp[3][gen]), accmean(&Fp[3][gen]), accmean(&SK2[3][gen]), accmean(&lociseg[3][gen]), accmean(&Qglobal[3][gen]), accmean(&meanQbottle[3][gen]), accmean(&HOglobal[3][gen]), accmean(&Hwglobal[3][gen]), accmean(&Hwbottle[3][gen]), accmean(&HT[3][gen]), accmean(&Fst[3][gen]), accmean(&Varq[3][gen]), accmean(&Varq_bt[3][gen]), accmean(&NeH[3][gen]), accmean(&NeH_bottle[3][gen]), accmean(&NeFst[3][gen]), accmean(&Nevq[3][gen]), accmean(&Nevq_bottle[3][gen]));		
		for (gen=genBPb; gen<=endbe; gen++)	fprintf(fgenBP,"%-6d  %-15f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-15f  %-15f  %-10f  %-15f\n", gen, accmean(&N[5][gen]), accmean(&nbott[5][gen]), accmean(&Ne[5][gen]), accmean(&fp[5][gen]), accmean(&Fp[5][gen]), accmean(&SK2[5][gen]), accmean(&lociseg[5][gen]), accmean(&Qglobal[5][gen]), accmean(&meanQbottle[5][gen]), accmean(&HOglobal[5][gen]), accmean(&Hwglobal[5][gen]), accmean(&Hwbottle[5][gen]), accmean(&HT[5][gen]), accmean(&Fst[5][gen]), accmean(&Varq[5][gen]), accmean(&Varq_bt[5][gen]), accmean(&NeH[5][gen]), accmean(&NeH_bottle[5][gen]), accmean(&NeFst[5][gen]), accmean(&Nevq[5][gen]), accmean(&Nevq_bottle[5][gen]));
		fprintf(fgenBP,"\n");			
	}
	if ((experiment==1.2)||(experiment==99))
	{
		fprintf(fgenBP,"Bottleneck (expA)\n");
		fprintf(fgenBP,"%-6s  %-15s  %-10s  %-15s  %-10s  %-10s  %-10s  %-15s  %-10s  %-10s  %-10s  %-10s  %-15s  %-10s  %-10s  %-10s  %-15s  %-10s  %-15s  %-15s  %-10s  %-15s\n", "gen", "N", "bottles", "NeVk", "fp", "Fp", "SK2", "loci_seg", "Qglobal", "Qbottle", "HO", "Hw", "Hw(bottle)", "HT", "Fst", "Varq", "Varq(bottle)", "NeH", "NeH(bottle)", "NeFst", "NeVq", "NeVq(bottle)");
		for (gen=0; gen<init; gen++)		fprintf(fgenBP,"%-6d  %-15f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-15f  %-15f  %-10f  %-15f\n", gen, accmean(&N[0][gen]), accmean(&nbott[0][gen]), accmean(&Ne[0][gen]), accmean(&fp[0][gen]), accmean(&Fp[0][gen]), accmean(&SK2[0][gen]), accmean(&lociseg[0][gen]), accmean(&Qglobal[0][gen]), accmean(&meanQbottle[0][gen]), accmean(&HOglobal[0][gen]), accmean(&Hwglobal[0][gen]), accmean(&Hwbottle[0][gen]), accmean(&HT[0][gen]), accmean(&Fst[0][gen]), accmean(&Varq[0][gen]), accmean(&Varq_bt[0][gen]), accmean(&NeH[0][gen]), accmean(&NeH_bottle[0][gen]), accmean(&NeFst[0][gen]), accmean(&Nevq[0][gen]), accmean(&Nevq_bottle[0][gen]));
		for (gen=init; gen<genBPb; gen++)	fprintf(fgenBP,"%-6d  %-15f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-15f  %-15f  %-10f  %-15f\n", gen, accmean(&N[3][gen]), accmean(&nbott[3][gen]), accmean(&Ne[3][gen]), accmean(&fp[3][gen]), accmean(&Fp[3][gen]), accmean(&SK2[3][gen]), accmean(&lociseg[3][gen]), accmean(&Qglobal[3][gen]), accmean(&meanQbottle[3][gen]), accmean(&HOglobal[3][gen]), accmean(&Hwglobal[3][gen]), accmean(&Hwbottle[3][gen]), accmean(&HT[3][gen]), accmean(&Fst[3][gen]), accmean(&Varq[3][gen]), accmean(&Varq_bt[3][gen]), accmean(&NeH[3][gen]), accmean(&NeH_bottle[3][gen]), accmean(&NeFst[3][gen]), accmean(&Nevq[3][gen]), accmean(&Nevq_bottle[3][gen]));		
		for (gen=genBPb; gen<=endbm; gen++)	fprintf(fgenBP,"%-6d  %-15f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-15f  %-15f  %-10f  %-15f\n", gen, accmean(&N[4][gen]), accmean(&nbott[4][gen]), accmean(&Ne[4][gen]), accmean(&fp[4][gen]), accmean(&Fp[4][gen]), accmean(&SK2[4][gen]), accmean(&lociseg[4][gen]), accmean(&Qglobal[4][gen]), accmean(&meanQbottle[4][gen]), accmean(&HOglobal[4][gen]), accmean(&Hwglobal[4][gen]), accmean(&Hwbottle[4][gen]), accmean(&HT[4][gen]), accmean(&Fst[4][gen]), accmean(&Varq[4][gen]), accmean(&Varq_bt[4][gen]), accmean(&NeH[4][gen]), accmean(&NeH_bottle[4][gen]), accmean(&NeFst[4][gen]), accmean(&Nevq[4][gen]), accmean(&Nevq_bottle[4][gen]));
		fprintf(fgenBP,"\n");			
	}

	//EXPERIMENT B
	if (experiment==0.2)
	{
		fprintf(fgenL,"Base Population and formation of Line\n");
		fprintf(fgenL,"%-6s  %-15s  %-10s  %-15s  %-10s  %-10s  %-10s  %-15s  %-10s  %-10s  %-10s  %-10s  %-15s  %-10s  %-10s  %-10s  %-15s  %-10s  %-15s  %-15s  %-10s  %-15s\n", "gen", "N", "bottles", "NeVk", "fp", "Fp", "SK2", "loci_seg", "Qglobal", "Qbottle", "HO", "Hw", "Hw(bottle)", "HT", "Fst", "Varq", "Varq(bottle)", "NeH", "NeH(bottle)", "NeFst", "NeVq", "NeVq(bottle)");
		for (gen=0; gen<tLINES; gen++)		fprintf(fgenL,"%-6d  %-15f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-15f  %-15f  %-10f  %-15f\n", gen, accmean(&N[0][gen]), accmean(&nbott[0][gen]), accmean(&Ne[0][gen]), accmean(&fp[0][gen]), accmean(&Fp[0][gen]), accmean(&SK2[0][gen]), accmean(&lociseg[0][gen]), accmean(&Qglobal[0][gen]), accmean(&meanQbottle[0][gen]), accmean(&HOglobal[0][gen]), accmean(&Hwglobal[0][gen]), accmean(&Hwbottle[0][gen]), accmean(&HT[0][gen]), accmean(&Fst[0][gen]), accmean(&Varq[0][gen]), accmean(&Varq_bt[0][gen]), accmean(&NeH[0][gen]), accmean(&NeH_bottle[0][gen]), accmean(&NeFst[0][gen]), accmean(&Nevq[0][gen]), accmean(&Nevq_bottle[0][gen]));
		for (gen=tLINES; gen<init; gen++)	fprintf(fgenL,"%-6d  %-15f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-15f  %-15f  %-10f  %-15f\n", gen, accmean(&N[1][gen]), accmean(&nbott[1][gen]), accmean(&Ne[1][gen]), accmean(&fp[1][gen]), accmean(&Fp[1][gen]), accmean(&SK2[1][gen]), accmean(&lociseg[1][gen]), accmean(&Qglobal[1][gen]), accmean(&meanQbottle[1][gen]), accmean(&HOglobal[1][gen]), accmean(&Hwglobal[1][gen]), accmean(&Hwbottle[1][gen]), accmean(&HT[1][gen]), accmean(&Fst[1][gen]), accmean(&Varq[1][gen]), accmean(&Varq_bt[1][gen]), accmean(&NeH[1][gen]), accmean(&NeH_bottle[1][gen]), accmean(&NeFst[1][gen]), accmean(&Nevq[1][gen]), accmean(&Nevq_bottle[1][gen]));
		fprintf(fgenL,"\n");
	}
	if ((experiment==2)||(experiment==99))
	{
		fprintf(fgenL,"Line Expansion (expB)\n");
		fprintf(fgenL,"%-6s  %-15s  %-10s  %-15s  %-10s  %-10s  %-10s  %-15s  %-10s  %-10s  %-10s  %-10s  %-15s  %-10s  %-10s  %-10s  %-15s  %-10s  %-15s  %-15s  %-10s  %-15s\n", "gen", "N", "bottles", "NeVk", "fp", "Fp", "SK2", "loci_seg", "Qglobal", "Qbottle", "HO", "Hw", "Hw(bottle)", "HT", "Fst", "Varq", "Varq(bottle)", "NeH", "NeH(bottle)", "NeFst", "NeVq", "NeVq(bottle)");
		for (gen=0; gen<tLINES; gen++)		fprintf(fgenL,"%-6d  %-15f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-15f  %-15f  %-10f  %-15f\n", gen, accmean(&N[0][gen]), accmean(&nbott[0][gen]), accmean(&Ne[0][gen]), accmean(&fp[0][gen]), accmean(&Fp[0][gen]), accmean(&SK2[0][gen]), accmean(&lociseg[0][gen]), accmean(&Qglobal[0][gen]), accmean(&meanQbottle[0][gen]), accmean(&HOglobal[0][gen]), accmean(&Hwglobal[0][gen]), accmean(&Hwbottle[0][gen]), accmean(&HT[0][gen]), accmean(&Fst[0][gen]), accmean(&Varq[0][gen]), accmean(&Varq_bt[0][gen]), accmean(&NeH[0][gen]), accmean(&NeH_bottle[0][gen]), accmean(&NeFst[0][gen]), accmean(&Nevq[0][gen]), accmean(&Nevq_bottle[0][gen]));
		for (gen=tLINES; gen<init; gen++)	fprintf(fgenL,"%-6d  %-15f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-15f  %-15f  %-10f  %-15f\n", gen, accmean(&N[1][gen]), accmean(&nbott[1][gen]), accmean(&Ne[1][gen]), accmean(&fp[1][gen]), accmean(&Fp[1][gen]), accmean(&SK2[1][gen]), accmean(&lociseg[1][gen]), accmean(&Qglobal[1][gen]), accmean(&meanQbottle[1][gen]), accmean(&HOglobal[1][gen]), accmean(&Hwglobal[1][gen]), accmean(&Hwbottle[1][gen]), accmean(&HT[1][gen]), accmean(&Fst[1][gen]), accmean(&Varq[1][gen]), accmean(&Varq_bt[1][gen]), accmean(&NeH[1][gen]), accmean(&NeH_bottle[1][gen]), accmean(&NeFst[1][gen]), accmean(&Nevq[1][gen]), accmean(&Nevq_bottle[1][gen]));
		for (gen=init; gen<=endLe; gen++)	fprintf(fgenL,"%-6d  %-15f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-15f  %-15f  %-10f  %-15f\n", gen, accmean(&N[6][gen]), accmean(&nbott[6][gen]), accmean(&Ne[6][gen]), accmean(&fp[6][gen]), accmean(&Fp[6][gen]), accmean(&SK2[6][gen]), accmean(&lociseg[6][gen]), accmean(&Qglobal[6][gen]), accmean(&meanQbottle[6][gen]), accmean(&HOglobal[6][gen]), accmean(&Hwglobal[6][gen]), accmean(&Hwbottle[6][gen]), accmean(&HT[6][gen]), accmean(&Fst[6][gen]), accmean(&Varq[6][gen]), accmean(&Varq_bt[6][gen]), accmean(&NeH[6][gen]), accmean(&NeH_bottle[6][gen]), accmean(&NeFst[6][gen]), accmean(&Nevq[6][gen]), accmean(&Nevq_bottle[6][gen]));
		fprintf(fgenL,"\n");
	}

	//EXPERIMENT C
	if ((experiment==3)||(experiment==99))			
	{
		fprintf(fgenmL,"Crossing lines (expC)\n");
		fprintf(fgenmL,"**LINE 1\n");
		fprintf(fgenmL,"%-6s  %-15s  %-10s  %-15s  %-10s  %-10s  %-10s  %-15s  %-10s  %-10s  %-10s  %-10s  %-15s  %-10s  %-10s  %-10s  %-15s  %-10s  %-15s  %-15s  %-10s  %-15s\n", "gen", "N", "bottles", "NeVk", "fp", "Fp", "SK2", "loci_seg", "Qglobal", "HO", "Hw", "Hw(bottle)", "HT", "Fst", "Varq", "Varq(bottle)", "NeH", "NeH(bottle)", "NeFst", "NeVq", "NeVq(bottle)");
		for (gen=0; gen<tLINES; gen++)		fprintf(fgenmL,"%-6d  %-15f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-15f  %-15f  %-10f  %-15f\n", gen, accmean(&N[0][gen]), accmean(&nbott[0][gen]), accmean(&Ne[0][gen]), accmean(&fp[0][gen]), accmean(&Fp[0][gen]), accmean(&SK2[0][gen]), accmean(&lociseg[0][gen]), accmean(&Qglobal[0][gen]), accmean(&meanQbottle[0][gen]), accmean(&HOglobal[0][gen]), accmean(&Hwglobal[0][gen]), accmean(&Hwbottle[0][gen]), accmean(&HT[0][gen]), accmean(&Fst[0][gen]), accmean(&Varq[0][gen]), accmean(&Varq_bt[0][gen]), accmean(&NeH[0][gen]), accmean(&NeH_bottle[0][gen]), accmean(&NeFst[0][gen]), accmean(&Nevq[0][gen]), accmean(&Nevq_bottle[0][gen]));
		for (gen=tLINES; gen<init; gen++)	fprintf(fgenmL,"%-6d  %-15f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-15f  %-15f  %-10f  %-15f\n", gen, accmean(&N[1][gen]), accmean(&nbott[1][gen]), accmean(&Ne[1][gen]), accmean(&fp[1][gen]), accmean(&Fp[1][gen]), accmean(&SK2[1][gen]), accmean(&lociseg[1][gen]), accmean(&Qglobal[1][gen]), accmean(&meanQbottle[1][gen]), accmean(&HOglobal[1][gen]), accmean(&Hwglobal[1][gen]), accmean(&Hwbottle[1][gen]), accmean(&HT[1][gen]), accmean(&Fst[1][gen]), accmean(&Varq[1][gen]), accmean(&Varq_bt[1][gen]), accmean(&NeH[1][gen]), accmean(&NeH_bottle[1][gen]), accmean(&NeFst[1][gen]), accmean(&Nevq[1][gen]), accmean(&Nevq_bottle[1][gen]));

		fprintf(fgenmL,"\n**LINE 2\n");
		fprintf(fgenmL,"%-6s  %-15s  %-10s  %-15s  %-10s  %-10s  %-10s  %-15s  %-10s  %-10s  %-10s  %-10s  %-15s  %-10s  %-10s  %-10s  %-15s  %-10s  %-15s  %-15s  %-10s  %-15s\n", "gen", "N", "bottles", "NeVk", "fp", "Fp", "SK2", "loci_seg", "Qglobal", "HO", "Hw", "Hw(bottle)", "HT", "Fst", "Varq", "Varq(bottle)", "NeH", "NeH(bottle)", "NeFst", "NeVq", "NeVq(bottle)");
		for (gen=0; gen<tLINES; gen++)		fprintf(fgenmL,"%-6d  %-15f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-15f  %-15f  %-10f  %-15f\n", gen, accmean(&N[0][gen]), accmean(&nbott[0][gen]), accmean(&Ne[0][gen]), accmean(&fp[0][gen]), accmean(&Fp[0][gen]), accmean(&SK2[0][gen]), accmean(&lociseg[0][gen]), accmean(&Qglobal[0][gen]), accmean(&meanQbottle[0][gen]), accmean(&HOglobal[0][gen]), accmean(&Hwglobal[0][gen]), accmean(&Hwbottle[0][gen]), accmean(&HT[0][gen]), accmean(&Fst[0][gen]), accmean(&Varq[0][gen]), accmean(&Varq_bt[0][gen]), accmean(&NeH[0][gen]), accmean(&NeH_bottle[0][gen]), accmean(&NeFst[0][gen]), accmean(&Nevq[0][gen]), accmean(&Nevq_bottle[0][gen]));
		for (gen=tLINES; gen<init; gen++)	fprintf(fgenmL,"%-6d  %-15f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-15f  %-15f  %-10f  %-15f\n", gen, accmean(&N[2][gen]), accmean(&nbott[2][gen]), accmean(&Ne[2][gen]), accmean(&fp[2][gen]), accmean(&Fp[2][gen]), accmean(&SK2[2][gen]), accmean(&lociseg[2][gen]), accmean(&Qglobal[2][gen]), accmean(&meanQbottle[2][gen]), accmean(&HOglobal[2][gen]), accmean(&Hwglobal[2][gen]), accmean(&Hwbottle[2][gen]), accmean(&HT[2][gen]), accmean(&Fst[2][gen]), accmean(&Varq[2][gen]), accmean(&Varq_bt[2][gen]), accmean(&NeH[2][gen]), accmean(&NeH_bottle[2][gen]), accmean(&NeFst[2][gen]), accmean(&Nevq[2][gen]), accmean(&Nevq_bottle[2][gen]));

		fprintf(fgenmL,"\n**CROSSED LINE\n");
		fprintf(fgenmL,"%-6s  %-15s  %-10s  %-15s  %-10s  %-10s  %-10s  %-15s  %-10s  %-10s  %-10s  %-10s  %-15s  %-10s  %-10s  %-10s  %-15s  %-10s  %-15s  %-15s  %-10s  %-15s\n", "gen", "N", "bottles", "NeVk", "fp", "Fp", "SK2", "loci_seg", "Qglobal", "HO", "Hw", "Hw(bottle)", "HT", "Fst", "Varq", "Varq(bottle)", "NeH", "NeH(bottle)", "NeFst", "NeVq", "NeVq(bottle)");
		for (gen=init; gen<=endmL; gen++)	fprintf(fgenmL,"%-6d  %-15f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-15f  %-15f  %-10f  %-15f\n", gen, accmean(&N[7][gen]), accmean(&nbott[7][gen]), accmean(&Ne[7][gen]), accmean(&fp[7][gen]), accmean(&Fp[7][gen]), accmean(&SK2[7][gen]), accmean(&lociseg[7][gen]), accmean(&Qglobal[7][gen]), accmean(&meanQbottle[7][gen]), accmean(&HOglobal[7][gen]), accmean(&Hwglobal[7][gen]), accmean(&Hwbottle[7][gen]), accmean(&HT[7][gen]), accmean(&Fst[7][gen]), accmean(&Varq[7][gen]), accmean(&Varq_bt[7][gen]), accmean(&NeH[7][gen]), accmean(&NeH_bottle[7][gen]), accmean(&NeFst[7][gen]), accmean(&Nevq[7][gen]), accmean(&Nevq_bottle[7][gen]));
		fprintf(fgenmL,"\n");
	}

	//Fst between lines
	if (experiment==99)
	{		
		fprintf(fgenLines,"** Natural Population N=%d\n", NINDNP);
		fprintf(fgenLines,"** Base Population N=%d (%d bottles)\n",ib*bottBP,bottBP);
		fprintf(fgenLines,"** Formation of lines: tLINES = %d N = %d (%d bottles)\n",tLINES,ib*bottL,bottL); 
		fprintf(fgenLines,"** Mating system %d\n",matingsystem);

		fprintf(fgenLines,"%-6s  %-15s  %-10s  %-15s  %-10s  %-10s  %-10s  %-15s  %-10s  %-10s  %-10s  %-15s  %-10s  %-10s  %-10s\n", "gen", "N", "bottles", "NeVk", "fp", "Fp", "SK2", "loci_seg", "Qglobal", "HO", "Hw", "Hw(bottle)", "HT", "Fst", "NeFst");
		for (gen=0; gen<tLINES; gen++)		fprintf(fgenLines,"%-6d  %-15f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f\n", gen, accmean(&N[0][gen]), accmean(&nbott[0][gen]), accmean(&Ne[0][gen]), accmean(&fp[0][gen]), accmean(&Fp[0][gen]), accmean(&SK2[0][gen]), accmean(&lociseg[0][gen]), accmean(&Qglobal[0][gen]), accmean(&HOglobal[0][gen]), accmean(&Hwglobal[0][gen]), accmean(&Hwbottle[0][gen]), accmean(&HT[0][gen]), accmean(&Fst[0][gen]), accmean(&NeFst[0][gen]));
		for (gen=tLINES; gen<init; gen++)	fprintf(fgenLines,"%-6d  %-15f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f  %-15f  %-10f  %-10f  %-10f\n", gen, (accmean(&N[1][gen]) + accmean(&N[2][gen]))/(double)lines, (accmean(&nbott[1][gen]) + accmean(&nbott[2][gen]))/(double)lines, (accmean(&Ne[1][gen]) + accmean(&Ne[2][gen]))/(double)lines, (accmean(&fp[1][gen]) + accmean(&fp[2][gen]))/(double)lines, (accmean(&Fp[1][gen]) + accmean(&Fp[2][gen]))/(double)lines, (accmean(&SK2[1][gen]) + accmean(&SK2[2][gen]))/(double)lines, (accmean(&lociseg[1][gen]) + accmean(&lociseg[2][gen]))/(double)lines, (accmean(&Qglobal[1][gen]) + accmean(&Qglobal[2][gen]))/(double)lines, (accmean(&HOglobal[1][gen]) + accmean(&HOglobal[2][gen]))/(double)lines, (accmean(&Hwglobal[1][gen]) + accmean(&Hwglobal[2][gen]))/(double)lines, (accmean(&Hwbottle[1][gen]) + accmean(&Hwbottle[2][gen]))/(double)lines, accmean(&HTlines[gen]), accmean(&Fst_lines[gen]), accmean(&NeFst_lines[gen]));
	}
}

/* ********************************************************************* */

sum_outline()
{                                                                         
	fprintf(fsum,"\n                 >>START EXP\n");
	fprintf(fsum,"	         .                 ________ (1.1)\n");
	fprintf(fsum,"                 .                /        *\n");
	fprintf(fsum,"                 .               /\n");
	fprintf(fsum,"_________________.              /__________ (1.2)\n");
	fprintf(fsum," Base.P :  (0.1)*|             /           *\n");
	fprintf(fsum,"        :        |      (1.0) /\n");
	fprintf(fsum,"        :        |__________*/\n");
	fprintf(fsum,"        :        .\n");
	fprintf(fsum,"        :        .   _________  (2)\n");
	fprintf(fsum,"        :        .  /         *\n");
	fprintf(fsum,"        :        . /\n");
	fprintf(fsum,"        :________./                           0. None Experiment (mant.)\n");
	fprintf(fsum," LINES  :  (0.2)*.                            1. Bottleneck-Expansion Experiment (A)\n");
	fprintf(fsum,"        :        .                            2. Lines-Expansion Experiment (B)\n");
	fprintf(fsum,"        :_______ .                            3. Crossing-Lines Experiment (C)\n");
	fprintf(fsum,"        :  L1   |.____________  (3)\n");
	fprintf(fsum,"        :_______|.  MIXED L   *               * outputs: data.ped and data.map\n");
	fprintf(fsum,"           L2    .\n");                                                                       
	fprintf(fsum,"\n\n");
	fprintf(fsum,"\n\nSeed: %d\nExperiment: mant(0.1; 0.2), A(1.0; 1.1; 1.2), B(2), C(3), ALL(99): %f\n Mating system (0 circular, 1 offspring pool, 2 complete panmixia)): %d\nMigration between bottles: %f\nNIND per Bottle (max 80): %d\nNumber Bottles - Base Population (max 100): %d\nGenerations before experiments: %d\nFormation of lines (generation): %d\nNumber Bottles - Lines (max 100): %d\n(expA)   Bottles bottleneck - BP (max 100): %d\n(expA)   Iniciation recovery - BP (generation): %d\n(expA)   Stop recovery - BP (generation): %d\n(expA 1.2) Mantain size until generation: %d\n(expA 1.1) Bottles expansion (max 100): %d\n(expA 1.1) Expansion until generation: %d\n(expA 1.1) Mantain size until generation: %d\n(expB)   Bottles expansion (max 100): %d\n(expB)   Expansion until generation: %d\n(expB)   Mantain size until generation: %d\n(expC)   Mantain mixed line until generation: %d\nLength of genome in Morgans(99:FreeRecom): %f\nNumber of Chromosomes (min 1, max 50): %d\nChromosome X (no 0, yes 1): %d\nNCRO (max 2000): %d\nNLOCI (2-30): %d\nLambda (neutral model): %f\nSample size for data_files :%d\nReplicates: %d\n\n", seed, experiment, matingsystem, mig, ib, bottBP, init, tLINES, bottL, bottr, genBPr, genBPb, endbm, botte, genBPe, endbm, bottLe, genLe, endLe, endmL, L, CRO, CROX, NCRO, NLOCI, Lambda, NINDs, replicates);   
}

/* ********************************************************************* */

distribution_out()
{
	if ((experiment!=0.2)&&(experiment!=2)&&(experiment!=3)) 	fprintf(fgenBP,"\n\nDISTRIBUTION q\n");
	if ((experiment==0.2)||(experiment==2)||(experiment==99))	fprintf(fgenL,"\n\nDISTRIBUTION q\n");
	if ((experiment==3)||(experiment==99))				fprintf(fgenmL,"\n\nDISTRIBUTION q\n");

	/* **************** EXPA **************** */   

	if (init%tgen==0)	init = init;	
 	else			init = init + tgen;

	if (genBPb%tgen==0)	genBPb = genBPb;	
 	else			genBPb = genBPb + tgen;

	if (experiment==0.1)
	{
		fprintf (fgenBP, "\n%12s  ", "Gen"); 
		for (i=0; i<(init/tgen); i++)	fprintf (fgenBP, "%8d  ", i*tgen);
		fprintf (fgenBP, "\n");

		// distribution of q values

		fprintf (fgenBP, "%12s  ", "q=0.0");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_00[0][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.0-0.1");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_00_01[0][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.1-0.2");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_01_02[0][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.2-0.3");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_02_03[0][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.3-0.4");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_03_04[0][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.4-0.5");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_04_05[0][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.5-0.6");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_05_06[0][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.6-0.7");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_06_07[0][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.7-0.8");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_07_08[0][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.8-0.9");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_08_09[0][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.9-1.0");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_09_10[0][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=1.0");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_10[0][i])/replicates);
		fprintf (fgenBP, "\n\n");
	}
	if (experiment==1.0)
	{
		fprintf (fgenBP, "\n%12s  ", "Gen"); 
		for (i=0; i<(genBPr/tgen); i++)	fprintf (fgenBP, "%8d  ", i*tgen);
		fprintf (fgenBP, "\n");

		// distribution of q values

		fprintf (fgenBP, "%12s  ", "q=0.0");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_00[0][i])/replicates);
		for (i=(init/tgen); i<(genBPr/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_00[3][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.0-0.1");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_00_01[0][i])/replicates);
		for (i=(init/tgen); i<(genBPr/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_00_01[3][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.1-0.2");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_01_02[0][i])/replicates);
		for (i=(init/tgen); i<(genBPr/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_01_02[3][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.2-0.3");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_02_03[0][i])/replicates);
		for (i=(init/tgen); i<(genBPr/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_02_03[3][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.3-0.4");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_03_04[0][i])/replicates);
		for (i=(init/tgen); i<(genBPr/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_03_04[3][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.4-0.5");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_04_05[0][i])/replicates);
		for (i=(init/tgen); i<(genBPr/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_04_05[3][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.5-0.6");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_05_06[0][i])/replicates);
		for (i=(init/tgen); i<(genBPr/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_05_06[3][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.6-0.7");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_06_07[0][i])/replicates);
		for (i=(init/tgen); i<(genBPr/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_06_07[3][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.7-0.8");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_07_08[0][i])/replicates);
		for (i=(init/tgen); i<(genBPr/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_07_08[3][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.8-0.9");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_08_09[0][i])/replicates);
		for (i=(init/tgen); i<(genBPr/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_08_09[3][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.9-1.0");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_09_10[0][i])/replicates);
		for (i=(init/tgen); i<(genBPr/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_09_10[3][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=1.0");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_10[0][i])/replicates);
		for (i=(init/tgen); i<(genBPr/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_10[3][i])/replicates);
		fprintf (fgenBP, "\n\n");
	}
	if ((experiment==1.1)||(experiment==99))
	{
		fprintf(fgenBP,"\n *** EXPA - Bottleneck & Expansion ***\n");
		fprintf (fgenBP, "\n%12s  ", "Gen"); 
		for (i=0; i<=(endbe/tgen); i++)	fprintf (fgenBP, "%8d  ", i*tgen);
		fprintf (fgenBP, "\n");

		// distribution of q values

		fprintf (fgenBP, "%12s  ", "q=0.0");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_00[0][i])/replicates);
		for (i=(init/tgen); i<(genBPb/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_00[3][i])/replicates);
		for (i=(genBPb/tgen); i<=(endbe/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_00[4][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.0-0.1");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_00_01[0][i])/replicates);
		for (i=(init/tgen); i<(genBPb/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_00_01[3][i])/replicates);
		for (i=(genBPb/tgen); i<=(endbe/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_00_01[4][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.1-0.2");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_01_02[0][i])/replicates);
		for (i=(init/tgen); i<(genBPb/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_01_02[3][i])/replicates);
		for (i=(genBPb/tgen); i<=(endbe/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_01_02[4][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.2-0.3");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_02_03[0][i])/replicates);
		for (i=(init/tgen); i<(genBPb/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_02_03[3][i])/replicates);
		for (i=(genBPb/tgen); i<=(endbe/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_02_03[4][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.3-0.4");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_03_04[0][i])/replicates);
		for (i=(init/tgen); i<(genBPb/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_03_04[3][i])/replicates);
		for (i=(genBPb/tgen); i<=(endbe/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_03_04[4][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.4-0.5");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_04_05[0][i])/replicates);
		for (i=(init/tgen); i<(genBPb/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_04_05[3][i])/replicates);
		for (i=(genBPb/tgen); i<=(endbe/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_04_05[4][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.5-0.6");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_05_06[0][i])/replicates);
		for (i=(init/tgen); i<(genBPb/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_05_06[3][i])/replicates);
		for (i=(genBPb/tgen); i<=(endbe/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_05_06[4][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.6-0.7");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_06_07[0][i])/replicates);
		for (i=(init/tgen); i<(genBPb/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_06_07[3][i])/replicates);
		for (i=(genBPb/tgen); i<=(endbe/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_06_07[4][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.7-0.8");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_07_08[0][i])/replicates);
		for (i=(init/tgen); i<(genBPb/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_07_08[3][i])/replicates);
		for (i=(genBPb/tgen); i<=(endbe/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_07_08[4][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.8-0.9");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_08_09[0][i])/replicates);
		for (i=(init/tgen); i<(genBPb/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_08_09[3][i])/replicates);
		for (i=(genBPb/tgen); i<=(endbe/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_08_09[4][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.9-1.0");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_09_10[0][i])/replicates);
		for (i=(init/tgen); i<(genBPb/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_09_10[3][i])/replicates);
		for (i=(genBPb/tgen); i<=(endbe/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_09_10[4][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=1.0");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_10[0][i])/replicates);
		for (i=(init/tgen); i<(genBPb/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_10[3][i])/replicates);
		for (i=(genBPb/tgen); i<=(endbe/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_10[4][i])/replicates);
		fprintf (fgenBP, "\n\n");
	}
	if ((experiment==1.2)||(experiment==99))
	{
		fprintf(fgenBP,"\n *** EXPA - Bottleneck ***\n");
		fprintf (fgenBP, "\n%12s  ", "Gen"); 
		for (i=0; i<=(endbm/tgen); i++)	fprintf (fgenBP, "%8d  ", i*tgen);
		fprintf (fgenBP, "\n");

		// distribution of q values

		fprintf (fgenBP, "%12s  ", "q=0.0");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_00[0][i])/replicates);
		for (i=(init/tgen); i<(genBPb/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_00[3][i])/replicates);
		for (i=(genBPb/tgen); i<=(endbm/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_00[5][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.0-0.1");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_00_01[0][i])/replicates);
		for (i=(init/tgen); i<(genBPb/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_00_01[3][i])/replicates);
		for (i=(genBPb/tgen); i<=(endbm/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_00_01[5][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.1-0.2");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_01_02[0][i])/replicates);
		for (i=(init/tgen); i<(genBPb/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_01_02[3][i])/replicates);
		for (i=(genBPb/tgen); i<=(endbm/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_01_02[5][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.2-0.3");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_02_03[0][i])/replicates);
		for (i=(init/tgen); i<(genBPb/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_02_03[3][i])/replicates);
		for (i=(genBPb/tgen); i<=(endbm/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_02_03[5][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.3-0.4");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_03_04[0][i])/replicates);
		for (i=(init/tgen); i<(genBPb/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_03_04[3][i])/replicates);
		for (i=(genBPb/tgen); i<=(endbm/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_03_04[5][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.4-0.5");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_04_05[0][i])/replicates);
		for (i=(init/tgen); i<(genBPb/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_04_05[3][i])/replicates);
		for (i=(genBPb/tgen); i<=(endbm/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_04_05[5][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.5-0.6");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_05_06[0][i])/replicates);
		for (i=(init/tgen); i<(genBPb/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_05_06[3][i])/replicates);
		for (i=(genBPb/tgen); i<=(endbm/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_05_06[5][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.6-0.7");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_06_07[0][i])/replicates);
		for (i=(init/tgen); i<(genBPb/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_06_07[3][i])/replicates);
		for (i=(genBPb/tgen); i<=(endbm/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_06_07[5][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.7-0.8");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_07_08[0][i])/replicates);
		for (i=(init/tgen); i<(genBPb/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_07_08[3][i])/replicates);
		for (i=(genBPb/tgen); i<=(endbm/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_07_08[5][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.8-0.9");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_08_09[0][i])/replicates);
		for (i=(init/tgen); i<(genBPb/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_08_09[3][i])/replicates);
		for (i=(genBPb/tgen); i<=(endbm/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_08_09[5][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=0.9-1.0");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_09_10[0][i])/replicates);
		for (i=(init/tgen); i<(genBPb/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_09_10[3][i])/replicates);
		for (i=(genBPb/tgen); i<=(endbm/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_09_10[5][i])/replicates);
		fprintf (fgenBP, "\n");
		fprintf (fgenBP, "%12s  ", "q=1.0");
		for (i=0; i<(init/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_10[0][i])/replicates);
		for (i=(init/tgen); i<(genBPb/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_10[3][i])/replicates);
		for (i=(genBPb/tgen); i<=(endbm/tgen); i++) fprintf (fgenBP, "%8.2f  ", accsum(&q_10[5][i])/replicates);
		fprintf (fgenBP, "\n\n");
	}

	/* **************** EXPB **************** */   

	if (experiment==0.2)
	{
		fprintf(fgenL,"\n *** LINE ***\n");
		fprintf (fgenL, "\n%12s  ", "Gen"); 
		for (i=0; i<(init/tgen); i++)	fprintf (fgenL, "%8d  ", i*tgen);
		fprintf (fgenL, "\n");

		// distribution of q values

		fprintf (fgenL, "%12s  ", "q=0.0");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_00[0][i])/replicates);
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_00[1][i])/replicates);
		fprintf (fgenL, "\n");
		fprintf (fgenL, "%12s  ", "q=0.0-0.1");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_00_01[0][i])/replicates);
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_00_01[1][i])/replicates);
		fprintf (fgenL, "\n");
		fprintf (fgenL, "%12s  ", "q=0.1-0.2");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_01_02[0][i])/replicates);
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_01_02[1][i])/replicates);
		fprintf (fgenL, "\n");
		fprintf (fgenL, "%12s  ", "q=0.2-0.3");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_02_03[0][i])/replicates);
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_02_03[1][i])/replicates);
		fprintf (fgenL, "\n");
		fprintf (fgenL, "%12s  ", "q=0.3-0.4");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_03_04[0][i])/replicates);
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_03_04[1][i])/replicates);
		fprintf (fgenL, "\n");
		fprintf (fgenL, "%12s  ", "q=0.4-0.5");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_04_05[0][i])/replicates);
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_04_05[1][i])/replicates);
		fprintf (fgenL, "\n");
		fprintf (fgenL, "%12s  ", "q=0.5-0.6");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_05_06[0][i])/replicates);
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_05_06[1][i])/replicates);
		fprintf (fgenL, "\n");
		fprintf (fgenL, "%12s  ", "q=0.6-0.7");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_06_07[0][i])/replicates);
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_06_07[1][i])/replicates);
		fprintf (fgenL, "\n");
		fprintf (fgenL, "%12s  ", "q=0.7-0.8");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_07_08[0][i])/replicates);
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_07_08[1][i])/replicates);
		fprintf (fgenL, "\n");
		fprintf (fgenL, "%12s  ", "q=0.8-0.9");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_08_09[0][i])/replicates);
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_08_09[1][i])/replicates);
		fprintf (fgenL, "\n");
		fprintf (fgenL, "%12s  ", "q=0.9-1.0");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_09_10[0][i])/replicates);
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_09_10[1][i])/replicates);
		fprintf (fgenL, "\n");
		fprintf (fgenL, "%12s  ", "q=1.0");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_10[0][i])/replicates);
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_10[1][i])/replicates);
		fprintf (fgenL, "\n\n");
	}
	if ((experiment==2)||(experiment==99))
	{
		fprintf(fgenL,"\n *** LINE - Expansion ***\n");
		fprintf (fgenL, "\n%12s  ", "Gen"); 
		for (i=0; i<=(endLe/tgen); i++)	fprintf (fgenL, "%8d  ", i*tgen);
		fprintf (fgenL, "\n");

		// distribution of q values

		fprintf (fgenL, "%12s  ", "q=0.0");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_00[0][i])/replicates);
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_00[1][i])/replicates);
		for (i=(init/tgen); i<=(endLe/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_00[6][i])/replicates);
		fprintf (fgenL, "\n");
		fprintf (fgenL, "%12s  ", "q=0.0-0.1");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_00_01[0][i])/replicates);
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_00_01[1][i])/replicates);
		for (i=(init/tgen); i<=(endLe/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_00_01[6][i])/replicates);
		fprintf (fgenL, "\n");
		fprintf (fgenL, "%12s  ", "q=0.1-0.2");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_01_02[0][i])/replicates);
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_01_02[1][i])/replicates);
		for (i=(init/tgen); i<=(endLe/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_01_02[6][i])/replicates);
		fprintf (fgenL, "\n");
		fprintf (fgenL, "%12s  ", "q=0.2-0.3");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_02_03[0][i])/replicates);
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_02_03[1][i])/replicates);
		for (i=(init/tgen); i<=(endLe/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_02_03[6][i])/replicates);
		fprintf (fgenL, "\n");
		fprintf (fgenL, "%12s  ", "q=0.3-0.4");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_03_04[0][i])/replicates);
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_03_04[1][i])/replicates);
		for (i=(init/tgen); i<=(endLe/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_03_04[6][i])/replicates);
		fprintf (fgenL, "\n");
		fprintf (fgenL, "%12s  ", "q=0.4-0.5");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_04_05[0][i])/replicates);
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_04_05[1][i])/replicates);
		for (i=(init/tgen); i<=(endLe/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_04_05[6][i])/replicates);
		fprintf (fgenL, "\n");
		fprintf (fgenL, "%12s  ", "q=0.5-0.6");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_05_06[0][i])/replicates);
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_05_06[1][i])/replicates);
		for (i=(init/tgen); i<=(endLe/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_05_06[6][i])/replicates);
		fprintf (fgenL, "\n");
		fprintf (fgenL, "%12s  ", "q=0.6-0.7");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_06_07[0][i])/replicates);
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_06_07[1][i])/replicates);
		for (i=(init/tgen); i<=(endLe/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_06_07[6][i])/replicates);
		fprintf (fgenL, "\n");
		fprintf (fgenL, "%12s  ", "q=0.7-0.8");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_07_08[0][i])/replicates);
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_07_08[1][i])/replicates);
		for (i=(init/tgen); i<=(endLe/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_07_08[6][i])/replicates);
		fprintf (fgenL, "\n");
		fprintf (fgenL, "%12s  ", "q=0.8-0.9");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_08_09[0][i])/replicates);
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_08_09[1][i])/replicates);
		for (i=(init/tgen); i<=(endLe/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_08_09[6][i])/replicates);
		fprintf (fgenL, "\n");
		fprintf (fgenL, "%12s  ", "q=0.9-1.0");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_09_10[0][i])/replicates);
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_09_10[1][i])/replicates);
		for (i=(init/tgen); i<=(endLe/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_09_10[6][i])/replicates);
		fprintf (fgenL, "\n");
		fprintf (fgenL, "%12s  ", "q=1.0");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_10[0][i])/replicates);
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_10[1][i])/replicates);
		for (i=(init/tgen); i<=(endLe/tgen); i++) fprintf (fgenL, "%8.2f  ", accsum(&q_10[6][i])/replicates);
		fprintf (fgenL, "\n\n");
	}

	/* **************** EXPC **************** */   

	if ((experiment==3)||(experiment==99))
	{
		fprintf(fgenmL,"\n *** BP ***\n");
		fprintf (fgenmL, "\n%12s  ", "Gen"); 
		for (i=0; i<(tLINES/tgen); i++)	fprintf (fgenmL, "%8d  ", i*tgen);
		fprintf (fgenmL, "\n");

		// distribution of q values

		fprintf (fgenmL, "%12s  ", "q=0.0");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_00[0][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.0-0.1");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_00_01[0][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.1-0.2");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_01_02[0][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.2-0.3");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_02_03[0][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.3-0.4");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_03_04[0][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.4-0.5");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_04_05[0][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.5-0.6");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_05_06[0][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.6-0.7");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_06_07[0][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.7-0.8");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_07_08[0][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.8-0.9");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_08_09[0][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.9-1.0");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_09_10[0][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=1.0");
		for (i=0; i<(tLINES/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_10[0][i])/replicates);
		fprintf (fgenmL, "\n\n");

		fprintf(fgenmL,"\n *** LINE 1 ***\n");
		fprintf (fgenmL, "\n%12s  ", "Gen"); 
		for (i=(tLINES/tgen); i<(init/tgen); i++)	fprintf (fgenmL, "%8d  ", i*tgen);
		fprintf (fgenmL, "\n");

		// distribution of q values

		fprintf (fgenmL, "%12s  ", "q=0.0");
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_00[1][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.0-0.1");
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_00_01[1][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.1-0.2");
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_01_02[1][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.2-0.3");
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_02_03[1][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.3-0.4");
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_03_04[1][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.4-0.5");
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_04_05[1][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.5-0.6");
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_05_06[1][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.6-0.7");
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_06_07[1][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.7-0.8");
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_07_08[1][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.8-0.9");
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_08_09[1][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.9-1.0");
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_09_10[1][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=1.0");
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_10[1][i])/replicates);
		fprintf (fgenmL, "\n\n");

		fprintf(fgenmL,"\n *** LINE 2 ***\n");
		fprintf (fgenmL, "\n%12s  ", "Gen"); 
		for (i=(tLINES/tgen); i<(init/tgen); i++)	fprintf (fgenmL, "%8d  ", i*tgen);
		fprintf (fgenmL, "\n");

		// distribution of q values

		fprintf (fgenmL, "%12s  ", "q=0.0");
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_00[2][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.0-0.1");
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_00_01[2][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.1-0.2");
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_01_02[2][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.2-0.3");
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_02_03[2][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.3-0.4");
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_03_04[2][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.4-0.5");
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_04_05[2][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.5-0.6");
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_05_06[2][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.6-0.7");
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_06_07[2][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.7-0.8");
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_07_08[2][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.8-0.9");
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_08_09[2][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.9-1.0");
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_09_10[2][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=1.0");
		for (i=(tLINES/tgen); i<(init/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_10[2][i])/replicates);
		fprintf (fgenmL, "\n\n");

		fprintf(fgenmL,"\n *** MIXED LINE ***\n");
		fprintf (fgenmL, "\n%12s  ", "Gen"); 
		for (i=(init/tgen); i<=(endmL/tgen); i++)	fprintf (fgenmL, "%8d  ", i*tgen);
		fprintf (fgenmL, "\n");

		// distribution of q values

		fprintf (fgenmL, "%12s  ", "q=0.0");
		for (i=(init/tgen); i<=(endmL/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_00[7][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.0-0.1");
		for (i=(init/tgen); i<=(endmL/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_00_01[7][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.1-0.2");
		for (i=(init/tgen); i<=(endmL/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_01_02[7][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.2-0.3");
		for (i=(init/tgen); i<=(endmL/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_02_03[7][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.3-0.4");
		for (i=(init/tgen); i<=(endmL/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_03_04[7][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.4-0.5");
		for (i=(init/tgen); i<=(endmL/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_04_05[7][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.5-0.6");
		for (i=(init/tgen); i<=(endmL/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_05_06[7][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.6-0.7");
		for (i=(init/tgen); i<=(endmL/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_06_07[7][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.7-0.8");
		for (i=(init/tgen); i<=(endmL/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_07_08[7][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.8-0.9");
		for (i=(init/tgen); i<=(endmL/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_08_09[7][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=0.9-1.0");
		for (i=(init/tgen); i<=(endmL/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_09_10[7][i])/replicates);
		fprintf (fgenmL, "\n");
		fprintf (fgenmL, "%12s  ", "q=1.0");
		for (i=(init/tgen); i<=(endmL/tgen); i++) fprintf (fgenmL, "%8.2f  ", accsum(&q_10[7][i])/replicates);
		fprintf (fgenmL, "\n\n");
	}
}
/* ********************************************************************* */


