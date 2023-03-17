
void RunAll(){

char fname[9][100];

sprintf(fname[0],"main_PMMA_70MeV.mac");
sprintf(fname[1],"main_PMMA_90MeV.mac");
sprintf(fname[2],"main_PMMA_100MeV.mac");

sprintf(fname[3],"main_PE_70MeV.mac");
sprintf(fname[4],"main_PE_90MeV.mac");
sprintf(fname[5],"main_PE_100MeV.mac");

sprintf(fname[6],"Zebra_70MeV.mac");
sprintf(fname[7],"Zebra_90MeV.mac");
sprintf(fname[8],"Zebra_100MeV.mac");

char Command[1612];

for(int n=0;n<9;n++) //Ncoll
{ 
	sprintf(Command, "Gate %s",fname[n]);//set output
	system(Command);//start move
}//RunALL fill
//*/
}
