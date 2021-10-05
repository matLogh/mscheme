// function return true if the entry is not suitable, e.g. two particles of same spin
// or the same combination of spins already present in the table
bool IsFound(int *vec, int nentries, int nelements, int **tab)
{
	for(int i=0; i<nentries; i++)
	{
		int fold = 0;
		for(int j=0; j<nelements; j++)
		{
			for(int k=0; k<nelements; k++)
			{
				if(vec[j] == vec[k] && k>=j+1) return true;
				if(vec[j] == tab[k][i] && i<nentries-1) 
				{
					fold++;
					break;
				}
			}
		}
		if(fold == nelements) return true;
	}

	return false;
}

// main function calculates J values from all combination of spins
// where nparticles is number of particles and spin is their 2*spin
Int_t mscheme(int spin = 7, int nparticles = 4)
{
	const int nfields = spin+1, nentries = TMath::Power(nfields,nparticles);
	int vSpin[nfields], index[nparticles];
	int **tabAll = new int *[nparticles+2];

	for(int i=0;i<nfields;i++)
		vSpin[i] = spin-2*i;		// array of all available spin values

	for(int i=0;i<nparticles;i++)
	{
		index[i] = 0;
		tabAll[i] = new int[nentries];	// 2D array of all available combinations
	}
	tabAll[nparticles] = new int[nentries];		// row for M values
	tabAll[nparticles+1] = new int[nentries];	// row for absolute sums - not used

	for(int i=0;i<nentries;i++)				// filling of tabAll
	{
		int sumM = 0, sumJ = 0;
		for(int j=0; j<nparticles; j++)
		{
			if( index[j]%nfields == 0 && j<=nparticles-1 && index[j]>0)
			{
				index[j] = 0;
				index[j+1]++;
			}
			sumM += vSpin[index[j]];
			sumJ += TMath::Abs(vSpin[index[j]]);
			tabAll[j][i] = vSpin[index[j]];
		}
		index[0]++;
		tabAll[nparticles][i] = sumM;
		tabAll[nparticles+1][i] = sumJ;
	}

	std::vector<int*> vMscheme;		// resulting m-scheme is polished tabAll
	for(int i=0;i<nentries;i++)
	{
		int *vec = new int[nparticles+2];
		for(int j=0; j<nparticles+2; j++)
		{
			vec[j] = tabAll[j][i];
		}
		if(IsFound(vec,i+1,nparticles,tabAll) || vec[nparticles]<0) continue;	// if non suitable entry - no filling mscheme
		vMscheme.push_back(vec);
	}

	TH1I *hist = new TH1I("hist_J","histogram of resulting M states",100,0,100);	// histogram of resulting M states
	for(int i=0;i<vMscheme.size();i++)												// histogram fill and screen output
	{		
		if(i==0) 
		{
			for(int j=0; j<nparticles; j++)
			{
				std::cout <<" m"<< j+1 << "  ";
				if(spin>10) std::cout <<" ";
			}
			std::cout << "   M" << std::endl;
		}
		for(int j=nparticles-1; j>=0; j--)
		{
			if(vMscheme[i][j]>0) std::cout << " ";
			std::cout << vMscheme[i][j] <<"/2 ";
		}
		if(nparticles%2==0)
		{
			std::cout <<"-> "<< vMscheme[i][nparticles]/2 << std::endl;
			hist->Fill(vMscheme[i][nparticles]/2);
		}
		else 
		{
			std::cout << "-> " << vMscheme[i][nparticles] <<"/2"<< std::endl;
			hist->Fill(vMscheme[i][nparticles]);
		}
	}

	std::cout << std::endl <<"Multiplicity\tx\tJ"<< std::endl;		// results in form of how many times and which J can be found
	int fold = 0;
	for(int i=hist->GetXaxis()->GetNbins(); i>0; i--)
	{
		int bin = hist->GetBinContent(i);
		if(bin>fold)
		{
			if(nparticles%2==0) std::cout <<"\t"<< bin-fold <<"\tx\t"<< i-1 << std::endl;
			else std::cout <<"\t"<< bin-fold <<"\tx\t"<< i-1 << "/2" << std::endl;
			fold = bin;
		}
	}

	hist->Draw("hist");

	return 0;
}
