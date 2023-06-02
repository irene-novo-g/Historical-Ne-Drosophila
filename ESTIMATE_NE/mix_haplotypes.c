// mix_haplotypes.c

#include "libhdr"
#include <time.h>
#include <stdlib.h>

#define NN 1000 // Maximum 1000 SNPs in fixed_alleles
#define MM 300000 // Maximum 300000 SNPs in data.map

char ch;
int x, s, r;
int numSNP;
char allele1, allele2, allele3, allele4, newallele1, newallele2;

FILE *find1, *find2, *fout, *fprueba;


main()
{
	srand(time(NULL));
	mix_haplotypes();
	return(0);
}

/* **************************************************************************** */

mix_haplotypes()
{
	// ********** read .ped **********

	find1 = fopen ("data1.ped","r");
	find2 = fopen ("data2.ped","r");
	fout = fopen ("mixed_pair.ped","w");
	fprueba = fopen ("prueba","w");

	fscanf(find1,"%d", &x);		//sex
	fscanf(find2,"%d", &x);		//sex
	fscanf(find1,"%d", &x);		//-9
	fscanf(find2,"%d", &x);		//-9
	fprintf(fout,"1 IND1 0 0 1 -9 ");

	while (!feof(find1))
	{
		s++;
		fprintf(fprueba,"s = %d\n", s);

		// read alleles from individual 1
		fscanf(find1,"%c", &ch); //space
		fscanf(find1,"%c", &ch); //1st allele
		allele1 = ch;
		fprintf(fprueba,"allele1 = %c\n", allele1);
		fscanf(find1,"%c", &ch); //space
		fscanf(find1,"%c", &ch); //2nd allele
		allele2 = ch;
		fprintf(fprueba,"allele2 = %c\n", allele2);
		if (allele1 == allele2)
		{
			newallele1 = allele1;
			fprintf(fprueba,"newallele1 = %c\n", newallele1);
		}
		else
		{
			r = rand() % 2;
			fprintf(fprueba,"r = %d\n", r);
			if (r == 0)
			{
				newallele1 = allele1;
			}
			else if (r == 1)
			{
				newallele1 = allele2;
			}
			fprintf(fprueba,"newallele1 = %c\n", newallele1);
		}

		// run out of SNPs: break
		if (allele1 == '\n')
		{
			break;
		}
	
		// read alleles from individual 2
		fscanf(find2,"%c", &ch); //space
		fscanf(find2,"%c", &ch); //1st allele
		allele3 = ch;
		fprintf(fprueba,"allele3 = %c\n", allele3);
		fscanf(find2,"%c", &ch); //space
		fscanf(find2,"%c", &ch); //2nd allele
		allele4 = ch;
		fprintf(fprueba,"allele4 = %c\n", allele4);
		if (allele3 == allele4)
		{
			newallele2 = allele3;
			fprintf(fprueba,"newallele2 = %c\n", newallele2);
		}
		else
		{
			r = rand() % 2;
			fprintf(fprueba,"r = %d\n", r);
			if (r == 0)
			{
				newallele2 = allele3;
			}
			else if (r == 1)
			{
				newallele2 = allele4;
			}
			fprintf(fprueba,"newallele2 = %c\n", newallele2);
		}

		fprintf(fout,"%c %c ", newallele1, newallele2);
	}
	
	fprintf(fout,"\n");

	fclose(find1);
	fclose(find2);
	fclose(fout);

	return(0);
}

