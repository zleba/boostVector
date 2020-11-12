
// Structure including properites of the track
struct track {
   double d0, z0, phi0, tanlambda;
};

struct event {
   track mu0, mu1;
   double t;       // time of event
   void toMicroM() { //from [cm] to [um]
      mu0.d0 *= 1e4;
      mu1.d0 *= 1e4;
      mu0.z0 *= 1e4;
      mu1.z0 *= 1e4;
   }
};

vector<event> getEvents(TString fName)
{
   TChain *tr = new TChain("variables");
   tr->Add(fName);

   if(!tr) {
      cout << "Something with reading" << endl;
      exit(0);
   }

   vector<event> events;


   double time, timeFrac;

   event evt;

   tr->SetBranchAddress("mu0_d0", &evt.mu0.d0);
   tr->SetBranchAddress("mu1_d0", &evt.mu1.d0);
   tr->SetBranchAddress("mu0_z0", &evt.mu0.z0);
   tr->SetBranchAddress("mu1_z0", &evt.mu1.z0);

   tr->SetBranchAddress("mu0_tanlambda", &evt.mu0.tanlambda);
   tr->SetBranchAddress("mu1_tanlambda", &evt.mu1.tanlambda);


   tr->SetBranchAddress("mu0_phi0", &evt.mu0.phi0);
   tr->SetBranchAddress("mu1_phi0", &evt.mu1.phi0);


   tr->SetBranchAddress("eventTimeSeconds", &time);
   tr->SetBranchAddress("eventTimeSecondsFractionRemainder", &timeFrac);

   for(int i = 0; i < tr->GetEntries(); ++i) {
      tr->GetEntry(i);
      evt.toMicroM();

      evt.t = time + timeFrac;


      events.push_back(evt);
   }

   //sort by time
   sort(events.begin(), events.end(), [](event e1, event e2) {return e1.t < e2.t;});

   return events;
}


// get the angles of vector orthogonal to the mu-mu plane
pair<double, double> getAngles(const event &e)
{
   TVector3 v0(cos(e.mu0.phi0), sin(e.mu0.phi0), e.mu0.tanlambda);
   TVector3 v1(cos(e.mu1.phi0), sin(e.mu1.phi0), e.mu1.tanlambda);
   TVector3 n = v0.Cross(v1);

   double phi = n.Phi();
   double angle = M_PI/2 - n.Theta();

   // the tan(angle) is in mrad
   return {phi, 1e3*tan(angle)};
}



// Fit xz and yz boost-angle
void fitBoost(const vector<event> &evts)
{
   TGraph *gr = new TGraph();

   for(auto e : evts) {
      double phi, tanLambda;
      tie(phi, tanLambda) = getAngles(e);
      gr->SetPoint(gr->GetN(), phi, tanLambda);
   }

   TCanvas *c = new TCanvas("canName", "");
   gr->Draw("ap");
   gr->GetXaxis()->SetTitle("#phi_{n} [rad]");
   gr->GetYaxis()->SetTitle("tan #lambda_{n} [mrad]");
   TF1 *f = new TF1("fun", "-[0]*cos(x) - [1]*sin(x)", -M_PI, M_PI);
   gr->Fit(f);

   gr->SetMinimum(-500);
   gr->SetMaximum(+500);
   f->Draw("same");

}



void plotFit()
{

   vector<event> evts = getEvents("ntuple1797.root");

   fitBoost(evts);

}
