/*
    dO3 = k2*O - k3*NO*O3 - J4*O3

y =      ( 0   1   2  3   4   5    6   7)
vector = (dNO2,dNO,dO3,dO,dOH,dHO2,dCO,dO1D)

p=           (0, 1, 2, 3,  4, 5, 6, 7, 8, 9, 10)
parameters = (J1,k2,k3,J4,k5,k6,k7,k8,k9,emis.no, H2O)

J1,J4,emis.no are time-dependent!

J1 <- J1a + 1e-2*photons // p[0]
J4 <- J4a + 1e-5*photons // p[3]
emis.no <- 1e6*photons // p[9]

*/

void calc_jacobian(double y[_numeq_], double jac[_numeq_][_numeq_], double p[_numpar_])
{

    //dNO2 = -J1*NO2 +k3*NO*O3 +k8*HO2*NO - k9*OH*NO2
	//jac[0][0] = -J1-k9*OH; // dNO2/dN02
	jac[0][0] = -p[0]-p[8]*y[4]; // dNO2/dN02
	//jac[0][1] = k3*O3+k8*HO2; //dN02/dNO
	jac[0][1] = p[2]*y[2]+p[7]*y[5]; //dN02/dNO
	//jac[0][2] = k3*NO; // dNO2/dO3
	jac[0][2] = p[2]*y[1]; // dNO2/dO3
	jac[0][3] = 0; // dN02/dO
	//jac[0][4] = -k9*NO2; //dNO2/dOH
	jac[0][4] = -p[8]*y[0]; //dNO2/dOH
	//jac[0][5] = k8*NO; //dNO2/dHO2
	jac[0][5] = p[7]*y[1]; //dNO2/dHO2
	jac[0][6] = 0; //dN02/CO
	jac[0][7] = 0; //DNO2/dO1D


    //dNO = J1*NO2 - k3*O3*NO - k8*HO2*NO + emis.no
	//jac[1][0] = J1; // dNO/dN02
	jac[1][0] = p[0]; // dNO/dN02
	//jac[1][1] = -k3*O3-k8*HO2; //dN0/dNO
	jac[1][1] = -p[2]*y[2]-p[7]*y[5]; //dN0/dNO
	//jac[1][2] = -k3*NO; // dNO/dO3
	jac[1][2] = -p[2]*y[1]; // dNO/dO3
	jac[1][3] = 0; // dN0/dO
	jac[1][4] = 0; //dNO/dOH
	//jac[1][5] = -k8*NO; //dNO/dHO2
	jac[1][5] = -p[7]*y[1]; //dNO/dHO2
	jac[1][6] = 0; //dN0/CO
	jac[1][7] = 0; //DNO/dO1D

    //dO3 = k2*O - k3*NO*O3 - J4*O3
	jac[2][0] = 0; // dO3/dN02
	//jac[2][1] = -k3*O3; // dO3/dNO
	jac[2][1] = -p[2]*y[2]; // dO3/dNO
	//jac[2][2] = -k3*NO-J4; // dO3/dO3
	jac[2][2] = -p[2]*y[1]-p[3]; // dO3/dO3
	//jac[2][3] = k2; // dO3/dO
	jac[2][3] = p[1]; // dO3/dO
	jac[2][4] = 0; // dO3/dOH
	jac[2][5] = 0; //dO3/dHO2
	jac[2][6] = 0; //dO3/CO
	jac[2][7] = 0; //DO3/dO1D

    //dO = J1*NO2 - k2*O + k5*O1D
	//jac[3][0] = J1; // dO/dN02
	jac[3][0] = p[0]; // dO/dN02
	jac[3][1] = 0; // dO/dNO
	jac[3][2] = 0; // dO/dO3
	//jac[3][3] = -k2; // dO/dO
	jac[3][3] = -p[1]; // dO/dO
	jac[3][4] = 0; // dO/dOH
	jac[3][5] = 0; // dO/dHO2
	jac[3][6] = 0; // dO/CO
	//jac[3][7] = k5; // DO/dO1D
	jac[3][7] = p[4]; // DO/dO1D

    //dOH = -k7*OH*CO + 2*k6*O1D*H2O - k9*OH*NO2
	//jac[4][0] = -k9*OH; // dOH/dN02
	jac[4][0] = -p[8]*y[4]; // dOH/dN02
	jac[4][1] = 0; // dOH/dNO
	jac[4][2] = 0; // dOH/dO3
	jac[4][3] = 0; // dOH/dO
	//jac[4][4] = -k7*CO-k9*NO2; // dOH/dOH
	jac[4][4] = -p[6]*y[6]-p[8]*y[0]; // dOH/dOH
	jac[4][5] = 0; //dOH/dHO2
	//jac[4][6] = -k7*OH; //dOH/CO
	jac[4][6] = -p[6]*y[4]; //dOH/CO
	//jac[4][7] = 2*k6*H2O; //DOH/dO1D
	jac[4][7] = 2*p[5]*p[10]; //DOH/dO1D

    //dHO2 = k7*OH*CO - k8*HO2*NO
	jac[5][0] = 0; // dHO2/dN02
	//jac[5][1] = -k8*HO2; // dHO2/dNO
	jac[5][1] = -p[7]*y[5]; // dHO2/dNO
	jac[5][2] = 0; // dHO2/dO3
	jac[5][3] = 0; // dHO2/dO
	//jac[5][4] = k7*CO; // dHO2/dOH
	jac[5][4] = p[6]*y[6]; // dHO2/dOH
	//jac[5][5] = -k8*NO; //dHO2/dHO2
	jac[5][5] = -p[7]*y[1]; //dHO2/dHO2
	//jac[5][6] = k7*OH; //dHO2/CO
	jac[5][6] = p[6]*y[4]; //dHO2/CO
	jac[5][7] = 0; //DHO2/dO1D

    // dCO = -k7*OH*CO
	jac[6][0] = 0; // dCO/dN02
	jac[6][1] = 0; // dCO/dNO
	jac[6][2] = 0; // dCO/dO3
	jac[6][3] = 0; // dCO/dO
	//jac[6][4] = -k7*CO; // dCO/dOH
	jac[6][4] = -p[6]*y[6]; // dCO/dOH
	jac[6][5] = 0; //dCO/dHO2
	//jac[6][6] = -k7*OH; // dCO/CO
	jac[6][6] = -p[6]*y[4]; // dCO/CO
	jac[6][7] = 0; //DCO/dO1D

    //dO1D = J4*O3 - k5*O1D - k6*O1D*H2O
	jac[7][0] = 0; // dO1D/dN02
	jac[7][1] = 0; // dO1D/dNO
	//jac[7][2] = J4; // dO1D/dO3
	jac[7][2] = p[3]; // dO1D/dO3
	jac[7][3] = 0; // dO1D/dO
	jac[7][4] = 0; // dO1D/dOH
	jac[7][5] = 0; //dO1D/dHO2
	jac[7][6] = 0; // dO1D/CO
	//jac[7][7] = -k5-k6*H20; //DCO/dO1D
	jac[7][7] = -p[4]-p[5]*p[10]; //DCO/dO1D

}

/*
y =      ( 0   1   2  3   4   5    6   7)
vector = (dNO2,dNO,dO3,dO,dOH,dHO2,dCO,dO1D)

p=           (0, 1, 2, 3,  4, 5, 6, 7, 8, 9, 10)
parameters = (J1,k2,k3,J4,k5,k6,k7,k8,k9,emis.no, H2O)
*/

void sode_system(double t, double y[_numeq_], double yout[_numeq_], double p[_numpar_])
{
    //dNO2 = -J1*NO2 +k3*NO*O3 +k8*HO2*NO - k9*OH*NO2
	yout[0] = -p[0]*y[0] + p[2]*y[1]*y[2]+p[7]*y[5]*y[1]-p[8]*y[4]*y[0];
    //dNO = J1*NO2 - k3*O3*NO - k8*HO2*NO + emis.no
    yout[1] = p[0]*y[0] - p[2]*y[2]*y[1]-p[7]*y[5]*y[1]+p[9];
    //dO3 = k2*O - k3*NO*O3 - J4*O3
    yout[2] = p[1]*y[3] - p[2]*y[1]*y[2] - p[3]*y[2];
    //dO = J1*NO2 - k2*O + k5*O1D
    yout[3] = p[0]*y[0] - p[1]*y[3] + p[4]*y[7];
    //dOH = -k7*OH*CO + 2*k6*O1D*H2O - k9*OH*NO2
    yout[4] = -p[6]*y[4]*y[6] + 2*p[5]*y[7]*p[10];
    //dHO2 = k7*OH*CO - k8*HO2*NO
    yout[5] = p[6]*y[4]*y[6] - p[7]*y[5]*y[1];
    // dCO = -k7*OH*CO
    yout[6] = -p[6]*y[4]*y[6];
    //dO1D = J4*O3 - k5*O1D - k6*O1D*H2O
    yout[7] = p[3]*y[2] - p[4]*y[7] - p[5]*y[7]*p[10];
}

