
void RunAll(){

char fname[3][100];

sprintf(fname[0],"ASPET_70MeV.mac");
sprintf(fname[1],"ASPET_90MeV.mac");
sprintf(fname[2],"ASPET_100MeV.mac");


char Command[1612];

for(int n=0;n<3;n++) //Ncoll
{ 
	sprintf(Command, "Gate %s",fname[n]);//set output
	system(Command);//start move
}//RunALL fill
//*/
}
