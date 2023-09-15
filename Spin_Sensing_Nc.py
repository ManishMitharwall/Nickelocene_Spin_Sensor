import numpy as np
np.seterr(all='ignore')
import matplotlib.pyplot as plt
import sys

def Gfunction(w, beta):
    G = w / (1 - np.exp(-beta * w))
    for l in range(len(G)):
        if np.isnan(G[l]):
            G[l] = 1/beta
    return G

def pauli_matrices(spin):
    """Generate the Pauli Spin Matrices for a given spin"""
    s = float(spin)
    n = int(2*s+1)
    sx=np.empty([n,n])
    sy=np.empty([n,n],dtype=complex)
    sz=np.empty([n,n])
    for a in range(1,n+1):
        for b in range(1,n+1):
            sx[a-1,b-1] = ( 0.5*((a==b+1) + (a+1==b))*np.sqrt((s+1)*(a+b-1)-a*b))
            sy[a-1,b-1] = 1j*(0.5*((a==b+1) - (a+1==b))*np.sqrt((s+1)*(a+b-1)-a*b))
            sz[a-1,b-1] = (s+1-a)*(a==b)
    return sx,sy,sz

def diagonalize(H):
    """Diagionalize the gien Hermitian hamiltonian """
    eigen_Values, eigen_Vectors = np.linalg.eigh(H)
    sort = np.argsort(eigen_Values)
    eigen_Values = np.array(eigen_Values)[sort]
    eigen_Vectors = eigen_Vectors[:,sort]
    return eigen_Values, eigen_Vectors

def S2_diagonalize(H,S):
    """Eigen value of S^2 operator from the Hamiltonian H"""
    eigenValues,eigenVectors=diagonalize(H)
    solv=np.linalg.solve(eigenVectors,(np.dot(S,eigenVectors)))
    svec=np.diag(solv)
    sval =  [0.5*(-1+np.sqrt(1+4*i)) for i in  svec]
    svec= [np.real(x) for x in svec]
    sval= [np.real(x) for x in sval]
    return  svec, sval

def Nc_Spin_matrix(M_spin, Nsites):
    Dim = int(2*M_spin+1)
    Sx1,Sy1,Sz1=pauli_matrices(1)
    for i in range(Nsites):
        if i==0:
            Nc_Sx = np.kron( Sx1, np.eye(Dim))
            Nc_Sy = np.kron( Sy1, np.eye(Dim))
            Nc_Sz = np.kron( Sz1, np.eye(Dim))
        else:
            Nc_Sx = np.kron( Nc_Sx, np.eye(Dim))
            Nc_Sy = np.kron( Nc_Sy, np.eye(Dim))
            Nc_Sz = np.kron( Nc_Sz, np.eye(Dim))            
    return Nc_Sx, Nc_Sy, Nc_Sz

def M_Spin_matrix(M_spin, Nsites):
    Dim = int(2*M_spin+1)
    Sxm,Sym,Szm=pauli_matrices(M_spin)
    M_Sx = np.zeros((3*Dim**(Nsites), 3*Dim**(Nsites), Nsites))
    M_Sy = np.zeros((3*Dim**(Nsites), 3*Dim**(Nsites), Nsites), dtype=complex)
    M_Sz = np.zeros((3*Dim**(Nsites), 3*Dim**(Nsites), Nsites))
    for i in range(Nsites):
        if i == 0:
            M_Sx_i = np.kron(np.eye(3), Sxm)
            M_Sy_i = np.kron(np.eye(3), Sym)
            M_Sz_i = np.kron(np.eye(3), Szm)   
        else:
            M_Sx_i = np.kron(np.eye(3), np.eye(Dim))
            M_Sy_i = np.kron(np.eye(3), np.eye(Dim))
            M_Sz_i = np.kron(np.eye(3), np.eye(Dim))           
        for k in range(1,Nsites):           
            if i==k:
                M_Sx_i = np.kron(M_Sx_i, Sxm)
                M_Sy_i = np.kron(M_Sy_i, Sym)
                M_Sz_i = np.kron(M_Sz_i, Szm)   
            else:
                M_Sx_i = np.kron(M_Sx_i, np.eye(Dim))
                M_Sy_i = np.kron(M_Sy_i, np.eye(Dim))
                M_Sz_i = np.kron(M_Sz_i, np.eye(Dim))     
        M_Sx[:,:,i] = M_Sx_i
        M_Sy[:,:,i] = M_Sy_i
        M_Sz[:,:,i] = M_Sz_i
    return M_Sx, M_Sy, M_Sz

def build_Hamil(S_Ni, S_M, D_aniso , M_spin, Nsites, B_field, J1, J2, ciclo = 0):
    Nc_Sx, Nc_Sy, Nc_Sz = S_Ni["Nc_Sx"], S_Ni["Nc_Sy"], S_Ni["Nc_Sz"]
    M_Sx, M_Sy, M_Sz = S_M["M_Sx"], S_M["M_Sy"], S_M["M_Sz"] 
    Dim = int(2*M_spin+1)
    H0 = np.zeros((3*Dim**(Nsites), 3*Dim**(Nsites)))
    H0 += D_aniso[0]*(M_Sx[:,:,0]@M_Sx[:,:,0]) + D_aniso[1]*np.real(M_Sy[:,:,0]@M_Sy[:,:,0]) + D_aniso[2]*(M_Sz[:,:,0]@M_Sz[:,:,0])
    H0 += B_field * M_Sz[:,:,0]
    H0 += 4*(Nc_Sz@Nc_Sz)
    H0 += B_field * Nc_Sz
    for site in range(Nsites-1):
        H0 += J1*( (M_Sx[:,:,site] @ M_Sx[:,:,site+1]) + np.real(M_Sy[:,:,site] @ M_Sy[:,:,site+1]) + (M_Sz[:,:,site] @ M_Sz[:,:,site+1]))
        H0 += D_aniso[0]*(M_Sx[:,:,site+1] @ M_Sx[:,:,site+1]) + D_aniso[1]*np.real(M_Sy[:,:,site+1] @ M_Sy[:,:,site+1]) + D_aniso[2]*(M_Sz[:,:,site+1] @ M_Sz[:,:,site+1])
        H0 += B_field * M_Sz[:,:,site+1]
        if site+2 < Nsites:
            H0 += J2*( (M_Sx[:,:,site] @ M_Sx[:,:,site+2]) + np.real(M_Sy[:,:,site] @ M_Sy[:,:,site+2]) + (M_Sz[:,:,site] @ M_Sz[:,:,site+2]))
    if ciclo == 1 and Nsites>=3:
        H0 += J1*( (M_Sx[:,:,0] @ M_Sx[:,:,Nsites-1]) + np.real(M_Sy[:,:,0] @ M_Sy[:,:,Nsites-1]) + (M_Sz[:,:,0] @ M_Sz[:,:,Nsites-1]))
    return H0

def boltz_pop(States_DOS, beta, E):
    Q = 0
    p = np.zeros((States_DOS))
    for n in range(States_DOS):
        p[n] = np.exp(-1*beta*E[n])
        Q += p[n]
    return p / Q

def cal_d2IdV2(States_DOS, E, C, V, p, Share, beta):
    d2IdV2_num_row = 0*np.array(V[2:])
    for n in range(States_DOS):
        for m in range(len(E)):
            Iwp = Gfunction((V-(E[m]-E[n])),beta) + Gfunction((V+(E[m]-E[n])),-beta)
            dIwp = np.diff(Iwp)/np.diff(V)
            d2Iwp = np.diff(dIwp)/np.diff(V[1:])
            d2IdV2_num_row =  d2IdV2_num_row +  p[n]* ( (C[:,n].T @ Share @ C[:,m])**2 + (C[:,m].T @ Share @ C[:,n])**2 ) * d2Iwp
    return d2IdV2_num_row

def normalize_data(data):
    vmin = np.min(data)
    vmax = np.max(data)
    # Normalize the data to the range [0, 1]
    normalized_data = (data - vmin) / (vmax - vmin)
    return normalized_data

def run_code(Jr, H0, S_Ni, S_M, Temp, w, States_DOS, coff_Nc, coff_M ):
    Nc_Sx, Nc_Sy, Nc_Sz = S_Ni["Nc_Sx"], S_Ni["Nc_Sy"], S_Ni["Nc_Sz"]
    M_Sx, M_Sy, M_Sz = S_M["M_Sx"], S_M["M_Sy"], S_M["M_Sz"] 
    d2IdV2 = []
    for ijr in Jr:
        H = H0 + ijr * (Nc_Sx@M_Sx[:,:,0] + np.real(Nc_Sy@M_Sy[:,:,0]) + Nc_Sz@M_Sz[:,:,0])
        E,C=diagonalize(H)
        kb = 0.086173324
        if (Temp < 0.2):
            beta = 1/(kb * 0.2)
        else:
            beta = 1/(kb * Temp)
        V = w.copy()
        Share =   coff_Nc*(Nc_Sx - 1j*Nc_Sy) + coff_M*(M_Sx[:,:,0] - 1j*M_Sy[:,:,0])
        p = boltz_pop(States_DOS, beta, E)
        der_cur = cal_d2IdV2(States_DOS, E, C, V, p, Share, beta)
        d2IdV2.append(der_cur)
    return np.real(d2IdV2)

def parameter_print(M_spin, D_aniso, Temp, B_field, J_Nc, Nsites, J1, J2, States_DOS, energy, coff_Nc, coff_M ):
    header_pr = f"""\n\n\n\t\t\t\t--------\t\t\t\t\t\n
Input Parameters are :
Spin of Metal = {M_spin}
D for metal is {D_aniso}  
Temperature is {Temp}
Magnetic field is {B_field}
Number of sites are {Nsites}
Exchange coupling between Nc and metal is {J_Nc} 
Coupling between nearest neighbour is {J1}
Coupling with second nearest neighbour is {J2}
Number of states taken for boltzman distribution {States_DOS}
Ploting range of energy is {energy}
Coupling of nickelocene spin to the metallic tip appex is {coff_Nc}
Coupling of surface spin to the underlying substrate is {coff_M}
\n\n\t\t\t\t--------\t\t\t\t\t\n\n
"""
    print(header_pr)

def read_parameters(filename):

    def convert_value(value):
        # Try to convert the value to a float
        try:
            return float(value)
        except ValueError:
            # If conversion to float fails, assume it's a list
            # Remove leading and trailing square brackets and split by commas
            if '[]' in value:
                value = value.strip('[]')
                value = [float(item.strip()) for item in value.split(',')]
                return value
            if '/' in value:
                parts = value.split('/')
                if len(parts) == 2:
                    try:
                        numerator = float(parts[0])
                        denominator = float(parts[1])
                        return numerator / denominator
                    except ValueError:
                        pass
            # If it's neither a float nor a fraction, assume it's a list
            # Remove leading and trailing square brackets and split by commas
            value = value.strip('[]')
            value = [float(item.strip()) for item in value.split(',')]
            return value
        
    parameters = {
        "m_spin": None,
        "d_aniso": None,
        "temp": None,
        "b_field": None,
        "j_nc": None,
        "nsites": None,
        "j1": None,
        "j2": None,
        "states_dos": None,
        "energy": None,
        "contrast": None,
        "coeff_nc": None,
        "coeff_m": None,
        "ciclo": None
    }

    with open(filename, 'r') as file:
        for line in file:
            # Remove leading and trailing whitespaces
            line = line.strip()

            # Skip empty lines and lines starting with #
            if not line or line.startswith('#'):
                continue

            # Split the line at the first '#' character to remove comments
            line = line.split('#', 1)[0].strip()

            # Split the line into key and value (assuming they are separated by '=')
            parts = line.split('=')
            if len(parts) == 2:
                key = parts[0].strip().lower()  # Convert to lowercase
                value = parts[1].strip()

                # Check if the key is one of the specified parameters
                if key in parameters:
                    # Convert the value using the convert_value function
                    parameters[key] = convert_value(value)

    return parameters

class Nc_sensing:
    def __init__(self) -> None:
        self.M_spin = None
        self.D_aniso = None
        self.Temp = 0.5
        self.B_field = 0.0
        self.J_Nc = 3
        self.Nsites = 1
        self.J1 = 6
        self.J2 = 0
        self.States_DOS = 3
        self.energy = 15
        self.contrast = 0.7
        self.coff_Nc = 3 
        self.coff_M  = 1
        self.ciclo = 0
    
    def Nc_spin_matrix(self):
        Nc_Sx, Nc_Sy, Nc_Sz = Nc_Spin_matrix(self.M_spin, self.Nsites)
        self.S_Ni = {"Nc_Sx" : Nc_Sx,
                     "Nc_Sy" : Nc_Sy,
                     "Nc_Sz" : Nc_Sz,
                     }
        return self.S_Ni
        
    def M_Spin_matrix(self):
        M_Sx, M_Sy, M_Sz = M_Spin_matrix(self.M_spin, self.Nsites)
        self.S_M = {"M_Sx" : M_Sx,
                     "M_Sy" : M_Sy,
                     "M_Sz" : M_Sz,
                     }        
        return self.S_M

    def bulid_Hamil(self):
        self.H0 = build_Hamil(self.S_Ni, self.S_M, self.D_aniso , self.M_spin, self.Nsites, self.B_field*0.1, self.J1, self.J2, ciclo=self.ciclo)

    def check_parameters(self):
        if len(self.D_aniso)<3 or (not isinstance(self.D_aniso, list)) : raise TypeError("D_aniso must be in [] line [0,0,5]")
        if not self.M_spin: raise ValueError("Please give M_spin")

    def parameter_print(self):
        parameter_print(self.M_spin, self.D_aniso, self.Temp, self.B_field, self.J_Nc, self.Nsites, self.J1, self.J2, self.States_DOS, self.energy, self.coff_Nc, self.coff_M )
    
    def kernel(self):
        self.check_parameters()
        self.parameter_print()
        self.Nc_spin_matrix()
        self.M_Spin_matrix()
        self.bulid_Hamil()
        w = np.arange(-self.energy,self.energy+0.1,0.1)
        Jr = np.linspace(0,self.J_Nc, 100)
        self.d2IdV2 = run_code(Jr, self.H0, self.S_Ni, self.S_M, self.Temp, w, self.States_DOS, coff_Nc=self.coff_Nc, coff_M=self.coff_M )
        if self.contrast:
            datap = normalize_data(self.d2IdV2)
            plt.title(f"Spin={self.M_spin} D={self.D_aniso} T={self.Temp} B={self.B_field}", size=12)
            plt.imshow(datap,extent=[min(w),max(w),np.min(Jr),np.max(Jr)],vmin=(self.contrast+0.01)/2, vmax=1-(self.contrast-0.01)/2 ,aspect='auto')
            plt.savefig(f"Plot_{self.M_spin}spin_{self.D_aniso}_{self.Temp}_{self.B_field}B.png",bbox_inches='tight',pad_inches=0.2,dpi=300)
        else:
            plt.title(f"Spin={self.M_spin} D={self.D_aniso} T={self.Temp} B={self.B_field}", size=12)
            plt.imshow(self.d2IdV2,extent=[min(w),max(w),np.min(Jr),np.max(Jr)],vmin=-20, vmax=20 ,aspect='auto')
            plt.savefig(f"Spin={self.M_spin}_D={self.D_aniso}_T={self.Temp}_B={self.B_field}.png",bbox_inches='tight',pad_inches=0.2,dpi=300)
        return self.d2IdV2


if __name__ == "__main__":
    my = Nc_sensing()
    if len(sys.argv) < 2:
        print("Give parameter file after python file");exit()
    param = read_parameters(sys.argv[1])
    if param["m_spin"]:
        my.M_spin = param["m_spin"]
    if param["d_aniso"]:
        my.D_aniso = param["d_aniso"]
    if param["temp"]:
        my.Temp = param["temp"]
    if param["b_field"]:
        my.B_field = param["b_field"]
    if param["j_nc"]:
        my.J_Nc = param["j_nc"]
    if param["nsites"]:
        my.Nsites = int(param["nsites"])
    if param["j1"]:
        my.J1 = param["j1"]
    if param["j2"]:
        my.J2 = param["j2"]
    if param["states_dos"]:
        my.States_DOS = int(param["states_dos"])
    if param["energy"]:
        my.energy = param["energy"]
    if param["contrast"]:
        my.contrast = param["contrast"]
    if param["coeff_nc"]:
        my.coff_Nc = param["coeff_nc"]
    if param["coeff_m"]:
        my.coff_M = param["coeff_m"]
    if param["ciclo"]:
        my.ciclo = param["ciclo"]
  
    my.kernel()
    plt.show()

 

