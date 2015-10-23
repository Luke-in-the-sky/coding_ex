{
	// This Macro analyzes the noise levels in (voltage) measurements performed under 
	// different conditions (temperature and electrical current).
	//
	// We have 11 txt-input files. We want to monitor if some values of applied current/temperature 
	// have an impact on the noise level in the measured output current (2nd column), the measured 
	// temperature (4th column) or the measured voltage (5th column)
	// Each input file has data stored in a tabular fashion, 
	// but not all colums are of float type. This means we will have to read the files
	// line by line, select what we need and save it to a new, float-only TNtuple object
	// for analysis and plotting.

  
	// Allocate list of file names
	ifstream ListOfFiles;
	const int NumOfFiles = 11;
	char* FileName[NumOfFiles] = {
		"Inj8-1Meas4-3(10K 1e-4A).txt",
		"Inj8-1Meas4-3(10K 1e-5A).txt",
		"Inj8-1Meas4-3(10K 1e-6A).txt",
		"Inj8-1Meas4-3(10K 1e-7A).txt",
		"Inj8-1Meas4-3(50K 1e-4A).txt",
		"Inj8-1Meas4-3(50K 1e-5A).txt",
		"Inj8-1Meas4-3(50K 1e-6A).txt",
		"Inj8-1Meas4-3(50K 1e-7A).txt",
		"Inj8-1Meas4-3(100K 1e-4A).txt",
		"Inj8-1Meas4-3(100K 1e-5A).txt",
		"Inj8-1Meas4-3(100K 1e-6A).txt"
		};
		
		
	// Allocate values of Set Temperature
	float Temp[NumOfFiles] = {
		10,
		10,
		10,
		10,
		50,
		50,
		50,
		50,
		100,
		100,
		100,    
	};


	// Allocate values of Set Current
	float Current[NumOfFiles] = {
		1e-4,
		1e-5,
		1e-6,
		1e-7,
		1e-4,
		1e-5,
		1e-6,
		1e-7,
		1e-4,
		1e-5,
		1e-6,
	};

	
	// Open a TCanvas for plotting
	TCanvas *can = new TCanvas("NoiseAnalysis","NoiseAnalysis", 120,120,1200, 600);
	can->Divide(2,2);
	int padindex=1;
	can->cd(padindex++);

	
	// Allocate objects for storing the relevant data
	TNtuple *result = new TNtuple("result", "result", "current:temp:VoltMean:VoltRMS");
	TNtuple *AllData = new TNtuple("all", "all", "current:temp:volt");
	char buffer[300];
	char buffer2[300];
	string sbuffer;
	float fbuffer;
	TH1F* htemp;

	
	// Loop on the input files, select relevant information and save it to the corresponding TNtuple
	bool GettingOut=0;
	for (int file=0; file< NumOfFiles; file++){
		TNtuple *temp = new TNtuple ("temp", "temp", "value");
		ListOfFiles.open(FileName[file]);
		cout << file << " - " << FileName[file] << " - " << Current[file] << " - " << Temp[file] << " - " << endl;
		
		
		// skip the header 
		bool HeaderFinished=0;
		for (int i=0; i<10; i++) {
			ListOfFiles >> sbuffer;
			if (sbuffer == "V") HeaderFinished=1;
			if (sbuffer == "V" && HeaderFinished) break;
		}
		getline(ListOfFiles, sbuffer);
		
		
		// Select what the data of interest (in this case, the entry in the 5th column)
		while (ListOfFiles) {
			ListOfFiles >> fbuffer;
			ListOfFiles >> fbuffer;
			ListOfFiles >> fbuffer;
			ListOfFiles >> fbuffer;
			ListOfFiles >> fbuffer;
			temp->Fill(fbuffer);
			AllData->Fill(Current[file], Temp[file], fbuffer);
			getline(ListOfFiles, sbuffer);
		}
		
		
		// Compute Average and RMS for entries in this input-file, save to relevant TNtuple.
		sprintf(buffer, "htemp-%i", file); cout << buffer << endl;
		sprintf(buffer2, "value >> %s", buffer); cout << buffer2 << endl;
		temp->Draw(buffer2, "", "goff");
		if (file==0) temp->Draw(buffer2, "", "");
		else temp->Draw(buffer2, "", "goff");
		htemp = (TH1F*)gDirectory->Get(buffer);
		result->Fill(Current[file], Temp[file], htemp->GetMean(), htemp->GetRMS());

		ListOfFiles.close();
		cout << endl;
	}


	// Plot it out
	result->SetMarkerStyle(21);
	padindex++;
	can->cd(padindex++);
	result->Draw("VoltMean/current:temp:-1*log(current)/log(10)", "", "colz");
	can->cd(padindex++);
	result->Draw("VoltRMS/VoltMean:temp:-1*log(current)/log(10)", "", "colz");

}