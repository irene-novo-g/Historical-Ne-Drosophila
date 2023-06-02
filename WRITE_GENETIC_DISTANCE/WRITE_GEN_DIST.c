// WRITE_GEN_DIST.c

#include "libhdr"

#define NN 2000000 // Maximum 2000000 SNPs in reference .map

int s, ss, x;
double w;
int numSNPref;
int chr_ref[NN], pos_ref[NN];
double gen_ref[NN];
int chr, pos;
double gen;
int found;

FILE *fref, *fin, *fout, *flost;

main()
{
	readreference();
	write_gen_dist();

	return(0);
}

/* **************************************************************************** */

readreference()
{
	// ********** read reference .map**********

	fref = fopen ("dataref.map","r");

	while (!feof(fref))
	{
		s ++;
		fscanf(fref,"%d", &x);		//chr
		chr_ref[s] = x;
		fscanf(fref,"%lf", &w);		//gen dist
		gen_ref[s] = w;
		fscanf(fref,"%d", &x);		//pos
		pos_ref[s] = x;
	}
	numSNPref = s - 1;

	fclose(fref);

	return(0);
}

/* **************************************************************************** */


write_gen_dist()
{
	// OUTPUT FILE
	fout = fopen ("out.map","w");
	// LOST SNPS FILE
	flost = fopen ("lost_list","w");

	// ********** read .map without gen distance and write .map with gen distance **********

	fin = fopen ("data_no_gen_dist.map","r");
	
	ss = 0;
	while (!feof(fin))
	{
		ss++;
		fscanf(fin,"%d", &x);		//chr
		chr = x;
		fscanf(fin,"%lf", &w);		//gen dist == 0
		fscanf(fin,"%d", &x);		//pos
		pos = x;

		found = 0;

		for (s=1; s<=numSNPref; s++)
		{
			if ((chr == chr_ref[s]) && (pos == pos_ref[s]))
			{
				gen = gen_ref[s];
				fprintf(fout,"%d SNP%d %f %d\n", chr, ss, gen, pos);
				found = 1;
				break;
			}
		}
		
		if (found == 0)
		{
			fprintf(fout,"%d SNP%d 0 %d\n", chr, ss, pos);
			fprintf(flost,"SNP%d\n", ss);
		}
	}

	fclose(fin);
	fclose(fout);
	fclose(flost);

	return(0);
}

