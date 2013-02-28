/* ------------------------------------------------------------------------ */

void GenContextRate(double* ST, double *fitQuant)
{
  int i;
  
  /* note that the table is filled in according to the ordering of 
     each codon using C=0,G=1,T=2,A=3 (summing of each nucleotide base 4) */
    Q[0][0]=Q[21][1]=0;
  Q[0][1]=Q[21][0]=0.49949243;     /* CCC->G = GGG->C =	     0.49949243	   0.0625047414 */
  Q[0][2]=Q[21][3]=2.73763617;     /* CCC->T = GGG->A =	     2.73763617	    0.143017678 */
  Q[0][3]=Q[21][2]=0.699398969;    /* CCC->A = GGG->T =	    0.699398969	   0.0563264239 */

  Q[1][0]=Q[5][1]=0;
  Q[1][1]=Q[5][0]=1.53321257;      /* CCG->G = CGG->C =	     1.53321257	    0.201204884 */
  Q[1][2]=Q[5][3]=27.9591204;      /* CCG->T = CGG->A =	     27.9591204	     1.76079219 */
  Q[1][3]=Q[5][2]=2.17616842;      /* CCG->A = CGG->T =	     2.17616842	    0.196229507 */

  Q[2][0]=Q[53][1]=0;
  Q[2][1]=Q[53][0]=0.496575844;    /* CCT->G = AGG->C =	    0.496575844	   0.0619832303 */
  Q[2][2]=Q[53][3]=2.33972251;     /* CCT->T = AGG->A =	     2.33972251	    0.121542271 */
  Q[2][3]=Q[53][2]=0.503430052;    /* CCT->A = AGG->T =	    0.503430052	   0.0409264478 */

  Q[3][0]=Q[37][1]=0;
  Q[3][1]=Q[37][0]=0.238857392;    /* CCA->G = TGG->C =	    0.238857392	   0.0296975073 */
  Q[3][2]=Q[37][3]=1.92505027;     /* CCA->T = TGG->A =	     1.92505027	    0.100607092 */
  Q[3][3]=Q[37][2]=0.535881624;    /* CCA->A = TGG->T =	    0.535881624	   0.0429739213 */

  Q[8][0]=Q[29][1]=1.18987981;     /* CTC->C = GAG->G =	     1.18987981	   0.0802485747 */
  Q[8][1]=Q[29][0]=0.273487539;    /* CTC->G = GAG->C =	    0.273487539	   0.0248541951 */
  Q[8][2]=Q[29][3]=0;
  Q[8][3]=Q[29][2]=0.286478042;    /* CTC->A = GAG->T =	    0.286478042	   0.0373643526 */

  Q[9][0]=Q[13][1]=1.39099009;     /* CTG->C = CAG->G =	     1.39099009	   0.0920613451 */
  Q[9][1]=Q[13][0]=0.361540482;    /* CTG->G = CAG->C =	    0.361540482	   0.0316738257 */
  Q[9][2]=Q[13][3]=0;
  Q[9][3]=Q[13][2]=0.299323283;    /* CTG->A = CAG->T =	    0.299323283	   0.0391137394 */

  Q[10][0]=Q[61][1]=1.02833941;    /* CTT->C = AAG->G =	     1.02833941	    0.158679646 */
  Q[10][1]=Q[61][0]=0.364167124;   /* CTT->G = AAG->C =	    0.364167124	   0.0315162334 */
  Q[10][2]=Q[61][3]=0;
  Q[10][3]=Q[61][2]=0.287679999;   /* CTT->A = AAG->T =	    0.287679999	   0.0373356322 */

  Q[11][0]=Q[45][1]=1.7167007;     /* CTA->C = TAG->G =	      1.7167007	    0.115705586 */
  Q[11][1]=Q[45][0]=0.399102329;   /* CTA->G = TAG->C =	    0.399102329	   0.0350924536 */
  Q[11][2]=Q[45][3]=0;
  Q[11][3]=Q[45][2]=0.385257976;   /* CTA->A = TAG->T =	    0.385257976	   0.0463699731 */

  Q[16][0]=Q[20][1]=0;
  Q[16][1]=Q[20][0]=0.516881906;   /* GCC->G = GGC->C =	    0.516881906	   0.0623716585 */
  Q[16][2]=Q[20][3]=2.32553943;    /* GCC->T = GGC->A =	     2.32553943	    0.123795573 */
  Q[16][3]=Q[20][2]=0.820251095;   /* GCC->A = GGC->T =	    0.820251095	   0.0647617296 */

  Q[17][0]=Q[4][1]=0;
  Q[17][1]=Q[4][0]=1.02263003;     /* GCG->G = CGC->C =	     1.02263003	    0.139613338 */
  Q[17][2]=Q[4][3]=20.6242194;     /* GCG->T = CGC->A =	     20.6242194	     1.30762809 */
  Q[17][3]=Q[4][2]=2.4149485;      /* GCG->A = CGC->T =	      2.4149485	    0.214690513 */

  Q[18][0]=Q[52][1]=0;
  Q[18][1]=Q[52][0]=0.482919104;   /* GCT->G = AGC->C =	    0.482919104	    0.059977258 */
  Q[18][2]=Q[52][3]=2.05210672;    /* GCT->T = AGC->A =	     2.05210672	    0.107997463 */
  Q[18][3]=Q[52][2]=0.668939654;   /* GCT->A = AGC->T =	    0.668939654	   0.0524826283 */

  Q[19][0]=Q[36][1]=0;
  Q[19][1]=Q[36][0]=0.441493676;   /* GCA->G = TGC->C =	    0.441493676	   0.0551262982 */
  Q[19][2]=Q[36][3]=1.88613884;    /* GCA->T = TGC->A =	     1.88613884	   0.0998119177 */
  Q[19][3]=Q[36][2]=0.970041474;   /* GCA->A = TGC->T =	    0.970041474	   0.0746090746 */

  Q[24][0]=Q[28][1]=1.69930077;    /* GTC->C = GAC->G =	     1.69930077	    0.115632459 */
  Q[24][1]=Q[28][0]=0.407646891;   /* GTC->G = GAC->C =	    0.407646891	   0.0600120603 */
  Q[24][2]=Q[28][3]=0;
  Q[24][3]=Q[28][2]=0.269561112;   /* GTC->A = GAC->T =	    0.269561112	   0.0355116703 */
  
  Q[25][0]=Q[12][1]=1.34469001;    /* GTG->C = CAC->G =	     1.34469001	   0.0908752944 */
  Q[25][1]=Q[12][0]=0.551831865;   /* GTG->G = CAC->C =	    0.551831865	   0.0476620137 */
  Q[25][2]=Q[12][3]=0;
  Q[25][3]=Q[12][2]=0.370433946;   /* GTG->A = CAC->T =	    0.370433946	   0.0481160216 */

  Q[26][0]=Q[60][1]=1.31763022;    /* GTT->C = AAC->G =	     1.31763022	     0.15601332 */
  Q[26][1]=Q[60][0]=0.372943188;   /* GTT->G = AAC->C =	    0.372943188	   0.0328060083 */
  Q[26][2]=Q[60][3]=0;
  Q[26][3]=Q[60][2]=0.256943759;   /* GTT->A = AAC->T =	    0.256943759	   0.0334233011 */

  Q[27][0]=Q[44][1]=2.11758864;    /* GTA->C = TAC->G =	     2.11758864	    0.141304942 */
  Q[27][1]=Q[44][0]=0.463852185;   /* GTA->G = TAC->C =	    0.463852185	   0.0407363263 */
  Q[27][2]=Q[44][3]=0;
  Q[27][3]=Q[44][2]=0.440305305;   /* GTA->A = TAC->T =	    0.440305305	    0.056939079 */

  Q[32][0]=Q[23][1]=0;
  Q[32][1]=Q[23][0]=0.58789432;    /* TCC->G = GGA->C =	     0.58789432	   0.0726919854 */
  Q[32][2]=Q[23][3]=2.3573874;     /* TCC->T = GGA->A =	      2.3573874	    0.122449821 */
  Q[32][3]=Q[23][2]=0.537510863;   /* TCC->A = GGA->T =	    0.537510863	   0.0430134513 */

  Q[33][0]=Q[7][1]=0;
  Q[33][1]=Q[7][0]=2.16188927;    /* TCG->G = CGA->C =	     2.16188927	    0.273005161 */
  Q[33][2]=Q[7][3]=28.0250364;    /* TCG->T = CGA->A =	     28.0250364	     1.74626231 */
  Q[33][3]=Q[7][2]=2.25988637;    /* TCG->A = CGA->T =	     2.25988637	    0.205555463 */

  Q[34][0]=Q[55][1]=0;
  Q[34][1]=Q[55][0]=0.700028129;   /* TCT->G = AGA->C =	    0.700028129	    0.085349789 */
  Q[34][2]=Q[55][3]=1.63105992;    /* TCT->T = AGA->A =	     1.63105992	   0.0850239197 */
  Q[34][3]=Q[55][2]=0.717607869;   /* TCT->A = AGA->T =	    0.717607869	    0.054844017 */

  Q[35][0]=Q[39][1]=0;
  Q[35][1]=Q[39][0]=0.455218434;   /* TCA->G = TGA->C =	    0.455218434	   0.0557958017 */
  Q[35][2]=Q[39][3]=1.65302071;    /* TCA->T = TGA->A =	     1.65302071	   0.0858997454 */
  Q[35][3]=Q[39][2]=0.579401651;   /* TCA->A = TGA->T =	    0.579401651	   0.0450482365 */

  Q[40][0]=Q[31][1]=0.853409308;   /* TTC->C = GAA->G =	    0.853409308	    0.100800939 */
  Q[40][1]=Q[31][0]=0.348563298;   /* TTC->G = GAA->C =	    0.348563298	   0.0506457826 */
  Q[40][2]=Q[31][3]=0;
  Q[40][3]=Q[31][2]=0.169248103;   /* TTC->A = GAA->T =	    0.169248103	   0.0221323993 */

  Q[41][0]=Q[15][1]=1.11714996;    /* TTG->C = CAA->G =	     1.11714996	   0.0766791028 */
  Q[41][1]=Q[15][0]=0.52769856;    /* TTG->G = CAA->C =	     0.52769856	   0.0456031237 */
  Q[41][2]=Q[15][3]=0;
  Q[41][3]=Q[15][2]=0.233963923;   /* TTG->A = CAA->T =	    0.233963923	   0.0280155205 */

  Q[42][0]=Q[63][1]=0.634389279;   /* TTT->C = AAA->G =	    0.634389279	   0.0751238191 */
  Q[42][1]=Q[63][0]=0.466417327;   /* TTT->G = AAA->C =	    0.466417327	    0.066736479 */
  Q[42][2]=Q[63][3]=0;
  Q[42][3]=Q[63][2]=0.245824523;   /* TTT->A = AAA->T =	    0.245824523	   0.0289239461 */

  Q[43][0]=Q[47][1]=1.21094061;    /* TTA->C = TAA->G =	     1.21094061	   0.0807523256 */
  Q[43][1]=Q[47][0]=0.565983138;   /* TTA->G = TAA->C =	    0.565983138	   0.0813378593 */
  Q[43][2]=Q[47][3]=0;
  Q[43][3]=Q[47][2]=0.467116545;   /* TTA->A = TAA->T =	    0.467116545	   0.0550414739 */

  Q[48][0]=Q[22][1]=0;
  Q[48][1]=Q[22][0]=0.604700869;   /* ACC->G = GGT->C =	    0.604700869	   0.0727373959 */
  Q[48][2]=Q[22][3]=2.52685924;    /* ACC->T = GGT->A =	     2.52685924	    0.133287544 */
  Q[48][3]=Q[22][2]=0.980611605;   /* ACC->A = GGT->T =	    0.980611605	   0.0763230465 */

  Q[49][0]=Q[6][1]=0;
  Q[49][1]=Q[6][0]=1.45676903;     /* ACG->G = CGT->C =	     1.45676903	    0.186359435 */
  Q[49][2]=Q[6][3]=34.5600009;     /* ACG->T = CGT->A =	     34.5600009	     2.13781585 */
  Q[49][3]=Q[6][2]=1.88818421;     /* ACG->A = CGT->T =	     1.88818421	    0.169686994 */

  Q[50][0]=Q[54][1]=0;
  Q[50][1]=Q[54][0]=0.562701186;   /* ACT->G = AGT->C =	    0.562701186	   0.0670710346 */
  Q[50][2]=Q[54][3]=2.35960961;    /* ACT->T = AGT->A =	     2.35960961	    0.122867731 */
  Q[50][3]=Q[54][2]=0.684482246;   /* ACT->A = AGT->T =	    0.684482246	   0.0535710106 */

  Q[51][0]=Q[38][1]=0;
  Q[51][1]=Q[38][0]=0.431405681;   /* ACA->G = TGT->C =	    0.431405681	   0.0518034257 */
  Q[51][2]=Q[38][3]=2.52414709;    /* ACA->T = TGT->A =	     2.52414709	    0.270259276 */
  Q[51][3]=Q[38][2]=0.740541187;   /* ACA->A = TGT->T =	    0.740541187	   0.0571986924 */

  Q[56][0]=Q[30][1]=1.7003273;     /* ATC->C = GAT->G =	      1.7003273	     0.11425539 */
  Q[56][1]=Q[30][0]=0.353084353;   /* ATC->G = GAT->C =	    0.353084353	   0.0311602624 */
  Q[56][2]=Q[30][3]=0;
  Q[56][3]=Q[30][2]=0.378385383;   /* ATC->A = GAT->T =	    0.378385383	   0.0451553466 */

  Q[57][0]=Q[14][1]=2.31157935;    /* ATG->C = CAT->G =	     2.31157935	    0.152649141 */
  Q[57][1]=Q[14][0]=0.478959857;   /* ATG->G = CAT->C =	    0.478959857	   0.0410458909 */
  Q[57][2]=Q[14][3]=0;
  Q[57][3]=Q[14][2]=0.446167626;   /* ATG->A = CAT->T =	    0.446167626	   0.0528829961 */

  Q[58][0]=Q[62][1]=1.7428248;     /* ATT->C = AAT->G =	      1.7428248	    0.115291805 */
  Q[58][1]=Q[62][0]=0.343104895;   /* ATT->G = AAT->C =	    0.343104895	   0.0293896518 */
  Q[58][2]=Q[62][3]=0;
  Q[58][3]=Q[62][2]=0.34053702;    /* ATT->A = AAT->T =	     0.34053702	   0.0400966766 */

  Q[59][0]=Q[46][1]=2.53979648;    /* ATA->C = TAT->G =	     2.53979648	     0.16816596 */
  Q[59][1]=Q[46][0]=0.40975183;    /* ATA->G = TAT->C =	     0.40975183	   0.0352254313 */
  Q[59][2]=Q[46][3]=0;
  Q[59][3]=Q[46][2]=0.481870464;   /* ATA->A = TAT->T =	    0.481870464	   0.0566958508 */

  /* stationary frequencies */
  ST[ 0] = 0.18417;  ST[ 4] = 0.02956;  ST[ 8] = 0.40883;  ST[12] = 0.37744; 
  ST[ 1] = 0.02748;  ST[ 5] = 0.02748;  ST[ 9] = 0.47252;  ST[13] = 0.47252; 
  ST[ 2] = 0.19046;  ST[ 6] = 0.02954;  ST[10] = 0.39947;  ST[14] = 0.38054; 
  ST[ 3] = 0.27310;  ST[ 7] = 0.02037;  ST[11] = 0.26949;  ST[15] = 0.43704; 
  ST[16] = 0.20056;  ST[20] = 0.20056;  ST[24] = 0.29944;  ST[28] = 0.29944; 
  ST[17] = 0.02956;  ST[21] = 0.18417;  ST[25] = 0.37744;  ST[29] = 0.40883; 
  ST[18] = 0.20754;  ST[22] = 0.16788;  ST[26] = 0.35549;  ST[30] = 0.26909; 
  ST[19] = 0.19101;  ST[23] = 0.16298;  ST[27] = 0.17462;  ST[31] = 0.47139; 
  ST[32] = 0.16298;  ST[36] = 0.19101;  ST[40] = 0.47139;  ST[44] = 0.17462; 
  ST[33] = 0.02037;  ST[37] = 0.27310;  ST[41] = 0.43704;  ST[45] = 0.26949; 
  ST[34] = 0.14516;  ST[38] = 0.23543;  ST[42] = 0.38807;  ST[46] = 0.23134; 
  ST[35] = 0.22160;  ST[39] = 0.22160;  ST[43] = 0.27840;  ST[47] = 0.27840; 
  ST[48] = 0.16788;  ST[52] = 0.20754;  ST[56] = 0.26909;  ST[60] = 0.35549; 
  ST[49] = 0.02954;  ST[53] = 0.19046;  ST[57] = 0.38054;  ST[61] = 0.39947; 
  ST[50] = 0.20331;  ST[54] = 0.20331;  ST[58] = 0.29669;  ST[62] = 0.29669; 
  ST[51] = 0.23543;  ST[55] = 0.14516;  ST[59] = 0.23134;  ST[63] = 0.38807; 

  if(fitQuant != NULL){
    double tmp[100] = {2.376344e-08, 1.193009e-06, 1.178988e-05, 5.989328e-05,
		       0.0002112965, 0.0005918941, 0.001414082, 0.003006856,
		       0.005849384, 0.01060785, 0.01817542, 0.02971528,
		       0.04670643, 0.07099246, 0.1048329, 0.1509573, 0.2126220,
		       0.2936690, 0.3985884, 0.5325819, 
		       0.7016296, 0.9125598, 1.173120, 1.492050, 1.879162, 
		       2.345416, 2.903003, 3.56543, 4.347609, 5.265943, 
		       6.338424, 7.58473, 9.026324, 10.68656, 12.59079, 
		       14.76648, 17.24334, 20.05342, 23.23128, 26.8141, 
		       30.84184, 35.3574, 40.40677, 46.03922, 52.30752, 
		       59.26808, 66.98122, 75.51142, 84.92756, 95.3032, 
		       106.7169, 119.2526, 133.0000, 148.0548, 164.5196, 
		       182.5038, 202.1248, 223.5083, 246.7889, 272.1114, 
		       299.6312, 329.5158, 361.9455, 397.1151, 435.2355, 
		       476.5351, 521.2622, 569.6871, 622.1048, 678.8384, 
		       740.2426, 806.7082, 878.6668, 956.5975, 1041.034, 
		       1132.572, 1231.882, 1339.721, 1456.947, 1584.541, 
		       1723.63, 1875.515, 2041.716, 2224.016, 2424.532, 
		       2645.801, 2890.896, 3163.594, 3468.601, 3811.888, 
		       4201.179, 4646.709, 5162.438, 5768.08, 6492.744, 
		       7381.969, 8512.812, 10031.41, 12270.46, 16311.35};
    
    for(i=0; i<100; i++){
      fitQuant[i] = tmp[i];
    }
  }
}

