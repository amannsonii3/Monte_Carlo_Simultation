#include<bits/stdc++.h>
#include<time.h>
using namespace std;
#define ll long long
#define OR ||
#define newl "\n"
#define rep(i,a,b) for(i=a;i<b;i++)
const ll INF = 1e9+7;

ll atoms = 700;
double box_dim = 10.0; // Length = Breadth = Height = 10 --> simulation space dimensions
ll moves = 100000; // Total number of iterations to do (counting the accepted moves only)
double T = 1;  
double kT = 1;
double sigma = 1; 
double epsilon = 1;
double cutoff = 3*sigma;

//random number between left and right --> left always = 0, right = 700(max.)
double random_number(double left, double right)
{
    double range = right-left;
    return (range)*(double(rand())/RAND_MAX) + left;   
    //RAND_MAX is the maximum random number that can be generated = 2147483647 
    //rand() will generate a random number
    //therefore double(rand())/RAND_MAX will generate a random value from (0-1)
    //returning a random number which will be multiple of the range as required
}

//minimum image convention
double min_img(double x)
{
    while(x<-5.0)
        x += 10.0;

    while(x>5.0)
        x -= 10.0;

    return x;
}

/*

The minimum image convention is a technique used in molecular simulations to reduce computational overhead by considering only the nearest periodic image of each atom. This is because in periodic simulations, the simulation box is replicated infinitely in all directions, creating periodic images of each atom. By using the minimum image convention, the simulation can effectively reduce the number of atoms that need to be considered for each interaction.

Minimum image convention yaane ki hum nearest image se interact karate hain... Agar kuch kaafi left mein hai to right mein nearest wala use karenge.. Vice versa

*/

//periodic boundary condition
double pbc(double x)
{
    while(x>10.0)
        x -= 10.0;

    while(x<0)
        x += 10.0;

    return x;
}

/* 

To simulate a bulk system in a finite box by imposing periodicity along the edges of the box. This means that when a particle crosses one face of the simulation box, it reappears on the opposite face, as if the box were replicated infinitely in all three dimensions. This creates an infinite, periodic array of boxes that represent the bulk system. PBC can be implemented using different boundary conditions, such as the minimum image convention

#DIFFERENCE BETWEEN MIC AND PBC

In summary, PBC is a technique to simulate a bulk system in a finite simulation box, while MIC is a technique to handle interactions between particles in a periodic system, by considering only the nearest periodic image of each particle.

*/

//initial energy calculation
double energy_calc(double box[][3])
{
    double energy = 0;
    ll i,j;
    for(i=0;i<atoms;i++)
    {
        for(j=i+1;j<atoms;j++)
        {
            double delx = box[j][0] - box[i][0]; 
            double dely = box[j][1] - box[i][1];
            double delz = box[j][2] - box[i][2];

            delx = min_img(delx);  
            dely = min_img(dely);
            delz = min_img(delz);
            // ensuring that if the point j wrt point i is displaced more than -5 to +5; the image of the point j is taken instead

            double r_sq = (delx*delx) + (dely*dely) + (delz*delz); //square of distance between two points = (atoms/molecules)

            if(r_sq>cutoff*cutoff) // only calculating for atoms within the cutoff radius
                continue;

            double sigma_sq = sigma*sigma;
            double temp = (sigma_sq/r_sq);
            double term2 = pow(temp,3);
            double term1 = term2*term2;
            double LJ = 4.0*epsilon*(term1-term2); //U(r) = 4ε[(σ/r)^12 - (σ/r)^6]
            energy += LJ;// calculating the LJ Potential Energy between the two points(from point 1 with point 2,3....n and then point 2 with point 1,3,4...n and so on)
        }
    }
    return energy;
}

//energy change after each displacement
double energy_change_calc(double box[][3],double old_x, double old_y, double old_z, double prev_energy,ll random_atom) // random_atom is the generated random_atom which was displaced and its new interactions with the previously fixed moecules gives the new LJ potential energy
{
    double old_interactions = 0.0;
    double new_interactions = 0.0;
    ll i;
    for(i=0;i<atoms;i++)
    {
        if(i!=random_atom)
        {
            double delx = box[random_atom][0] - box[i][0];
            double dely = box[random_atom][1] - box[i][1];
            double delz = box[random_atom][2] - box[i][2];

            delx = min_img(delx);
            dely = min_img(dely);
            delz = min_img(delz);

            double r_sq = (delx*delx) + (dely*dely) + (delz*delz);

            if(r_sq>cutoff*cutoff)
                continue;

            double sigma_sq = sigma*sigma;
            double temp = (sigma_sq/r_sq);
            double term2 = pow(temp,3);
            double term1 = term2*term2;
            double LJ = 4.0*epsilon*(term1-term2); 
            new_interactions += LJ;
        }
    }

    for(i=0;i<atoms;i++)
    {
        if(i!=random_atom)
        {
            double delx = old_x - box[i][0];
            double dely = old_y - box[i][1];
            double delz = old_z - box[i][2];

            delx = min_img(delx);
            dely = min_img(dely);
            delz = min_img(delz);

            double r_sq = (delx*delx) + (dely*dely) + (delz*delz);

            if(r_sq>cutoff*cutoff)
                continue;

            double sigma_sq = sigma*sigma;
            double temp = (sigma_sq/r_sq);
            double term2 = pow(temp,3);
            double term1 = term2*term2;
            double LJ = 4.0*epsilon*(term1-term2);
            old_interactions += LJ;
        }
    }

    return prev_energy - old_interactions + new_interactions;
}

int main() 
{
    ios::sync_with_stdio(0);
    cin.tie(0);
    srand(time(0));
    double box[atoms][3]; // a 2D array storing information of the points with atoms = atom number and array idicies 0=x,1=y,2=z
    ll i,j;
    vector<double>energy;
    vector<double> index {0,0,0};
    
    //starting with a fixed configuration
    for(i=0;i<atoms;i++)
    {
        for(j=0;j<3;j++)
        {
            box[i][j] = (int32_t)((index[j]+0.5)*(10/9)); // using the factor of 10/9 to evenly space out the distribution of the atoms in the box of L=10*10*10 and 700 atoms

        }
        index[0] = index[0]+1; // moving x-coordinate x to x+1 as we move from i to i+1 atom
        if(index[0] == 9) //shifting to y-coordinate space and moving y to y+1 as we move from i to i+1 atom
        {
          index[0] = 0;
          index[1] = index[1]+1;
          if(index[1]==9)// repeating the similar process for z to z+1
          {
            index[1]=0;
            index[2] = index[2]+1;
          }
        }
    }
    
    /*

    TO PRINT THE INITIAL CONFIGURATION OF THE 700 ATOMS
    
    for(int i = 0; i<atoms; i++){
        for(int j = 0; j<3; j++){
            cout << box[i][j] << " ";
        }
        cout << endl;
    }
    
    */

    //pushing the initial finite probability configuration
    energy.push_back(energy_calc(box)); 
    cout<<energy.back()<<newl; 
    
    ll accept = 0;    
    
    //begining with the simulation
    while(1)
    {
        //selecting a random atom
        ll random_atom = (ll)(random_number(0,atoms));

        //storing old coordinates for using, if rejected
        double old_x = box[random_atom][0];
        double old_y = box[random_atom][1];
        double old_z = box[random_atom][2];

        //giving random displacement
        box[random_atom][0] += random_number(0,1.0) - 0.5;
        box[random_atom][1] += random_number(0,1.0) - 0.5;
        box[random_atom][2] += random_number(0,1.0) - 0.5;

        //applying pbc to new coordinates since the random displacement can go beyond the box dimensions
        box[random_atom][0] = pbc(box[random_atom][0]);
        box[random_atom][1] = pbc(box[random_atom][1]);
        box[random_atom][2] = pbc(box[random_atom][2]);

        double new_energy = energy_change_calc(box,old_x,old_y,old_z,energy.back(),random_atom);

        //energy.back() -> last valid configuration's energy
        double energy_change = new_energy - energy.back();

        //energy change less than 0 -> finite probability
        //so we accept the move
        if(new_energy<=energy.back())
        {
            energy.push_back(new_energy);
            
            accept++;
            //if(accept%5==0)
                cout<<energy.back()<<newl;
        }
        else
        {
            double check = exp(-energy_change/kT);

            //calling a random number between 0 and 1
            double prob = random_number(0.0,1.0);

            //if random number <= probability term -> finite probability
            //accept the move 
            //To try not to fall into a local minima and stop the loop only when we obtain the global minima
            if(prob<=check)
            {
                energy.push_back(new_energy);
                
                accept++;
                //if(accept%5==0)
                    cout<<energy.back()<<newl;
            }

            //else reject the move and restore the old configuration
            else
            {
                box[random_atom][0] = old_x;
                box[random_atom][1] = old_y;
                box[random_atom][2] = old_z;
            }
        }
        if(accept==moves)
        break;
    }
    
    /*

    TO PRINT THE FINAL CONFIGURATION OF THE 700 ATOMS
    
    for(int i = 0; i<atoms; i++){
        for(int j = 0; j<3; j++){
            cout << box[i][j] << " ";
        }
        cout << endl;
    }
    
    */

    return 0;
}
