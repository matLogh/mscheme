// function return true if the entry is not suitable, e.g. two particles of same spin
// or the same combination of spins already present in the table
bool IsFound(int *vec, int nelements, std::vector<int*> tab)
{
	if(tab.size()==0)
	{
		for(int j=0; j<nelements; j++)
			for(int k=0; k<nelements; k++)
				if(vec[j] == vec[k] && k>=j+1) return true;
	}

	for(int i=0; i<tab.size(); i++)
	{
		int fold = 0;
		for(int j=0; j<nelements; j++)
		{
			for(int k=0; k<nelements; k++)
			{
				if(vec[j] == vec[k] && k>=j+1) return true;
				if(vec[j] == tab[i][k])
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
	std::vector<int*> vMscheme;		// resulting m-scheme

	for(int i=0;i<nfields;i++)
		vSpin[i] = spin-2*i;		// array of all available spin values

	for(int i=0;i<nparticles;i++)
	{
		index[i] = 0;
	}

	std::cout << setprecision(3);
	for(int i=0;i<nentries;i++)				// filling of mscheme table
	{
		int sumM = 0, sumJ = 0;
		int *vec = new int[nparticles+2];
		for(int j=0; j<nparticles; j++)
		{
			if( index[j]%nfields == 0 && j<=nparticles-1 && index[j]>0)
			{
				index[j] = 0;
				index[j+1]++;
			}
			sumM += vSpin[index[j]];
			sumJ += TMath::Abs(vSpin[index[j]]);
			vec[j] = vSpin[index[j]];
		}
		index[0]++;
		vec[nparticles] = sumM;
		vec[nparticles+1] = sumJ;
		if(i%200000==0) std::cout << "Processed " << setw(7) << (double)i/((double)nentries)*100 << setw(3) << "%\r" << std::flush;
		if(IsFound(vec,nparticles,vMscheme) || vec[nparticles]<0)
		{
			delete[] vec;										// clear RAM space
			continue;											// if non suitable entry - no filling mscheme
		}
		vMscheme.push_back(vec);
	}
	std::cout << std::endl;

	TH1I *hist = new TH1I("hist_J","histogram of resulting M states",100,0,100);	// histogram of resulting M states
	for(int i=0;i<vMscheme.size();i++)												// histogram fill and screen output
	{		
		if(i==0) 
		{
			for(int j=0; j<nparticles; j++)
			{
				if(j==0) std::cout << setw(7) <<"m1";
				else std::cout << setw(7) <<"m"<< j+1;
			}
			std::cout << setw(7) << "M" << std::endl;
		}
		for(int j=nparticles-1; j>=0; j--)
		{
			std::cout << setw(5) << vMscheme[i][j] <<"/2 ";
		}
		if(nparticles%2==0)
		{
			std::cout << setw(4) <<"-> "<< vMscheme[i][nparticles]/2 << std::endl;
			hist->Fill(vMscheme[i][nparticles]/2);
		}
		else 
		{
			std::cout << setw(4) << "-> " << vMscheme[i][nparticles] <<"/2"<< std::endl;
			hist->Fill(vMscheme[i][nparticles]);
		}
	}

	std::cout << std::endl <<"Multiplicity  x     J"<< std::endl;		// results in form of how many times and which J can be found
	int fold = 0;
	for(int i=hist->GetXaxis()->GetNbins(); i>0; i--)
	{
		int bin = hist->GetBinContent(i);
		if(bin>fold)
		{
			std::cout << setw(6);
			if(nparticles%2==0) std::cout << bin-fold << setw(9) <<"x"<< setw(7) << i-1 << std::endl;
			else std::cout << bin-fold << setw(9) <<"x"<< setw(5) << i-1 << "/2" << std::endl;
			fold = bin;
		}
	}

	hist->Draw("hist");

	return 0;
}
