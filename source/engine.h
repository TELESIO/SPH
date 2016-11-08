#ifndef _SPH_ENGINE
#define _SPH_ENGINE

#include <vector>
#include "particle.h"
#include <random>
#include <cmath>
#include "solid_surface.h"
#include<glm/gtx/rotate_vector.hpp>
#include<unordered_map>
typedef particle::VEC VEC;






double Wpoly6(const double radiusSquared, const double h) {

    static double coefficient = 315.0/(64.0*M_PI*pow(h,9));
    static double hSquared = pow(h,2);

    return coefficient * pow(hSquared-radiusSquared, 3);
}

void Wpoly6Gradient(VEC& diffPosition, double radiusSquared, VEC& gradient , const double h) {

    static double coefficient = -945.0/(32.0*M_PI*pow(h,9));
    static double hSquared = pow(h,2);

    gradient =  coefficient * pow(hSquared-radiusSquared, 2) * diffPosition ;


}

inline double Wpoly6Laplacian(const double radiusSquared, const double h) {

    static double coefficient = -945.0/(32.0*M_PI*pow(h,9));
    static double hSquared = pow(h,2);

    return coefficient * (hSquared-radiusSquared) * (3.0*hSquared - 7.0*radiusSquared);
}

inline void WspikyGradient(VEC& diffPosition, double radiusSquared, VEC& gradient, const double h) {  //

    static double coefficient = -45.0/(M_PI*pow(h,6));

    double radius = sqrt(radiusSquared);

    gradient = coefficient * pow(h-radius, 2) * diffPosition / radius;
}


inline double WviscosityLaplacian(double radiusSquared, const double h) {

    static double coefficient = 45.0/(M_PI*pow(h,6));

    double radius = sqrt(radiusSquared);

    return coefficient * (h - radius);
}




//----------------------END OF UTILITY FUNCTIONS





constexpr const static  double STIFFNESS  = 3.0; // Nm/kg is gas constant of water vapor
constexpr const static  double REST_DENSITY = 998.29; //kg/m^3 is rest density of water particle
constexpr const static  double VISCOSITY = 3.5; // Ns/m^2 or Pa*s viscosity of water

constexpr const static  double  SURFACE_TENSION = 0.0728; // N/m
constexpr const static  double  SURFACE_THRESHOLD = 7.065;

constexpr const static  double GRAVITY_ACCELERATION = 9.80665;


constexpr const static  double WALL_DAMPING =-0.9; // wall damping constant
constexpr const static  double WALL_K = 10000.0; // wall spring constant



//template<class P>
//unsigned long hash(const P& particle , double radius){
//    typedef
//}


inline void hash_combine(std::size_t& seed) { return; }
template <typename T, typename... Rest>
inline void hash_combine(std::size_t& seed, const T& v, Rest... rest) {
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
    hash_combine(seed, rest...);
}


namespace std {
template<typename a, typename b>
struct hash< std::pair<a, b> > {
private:
    const hash<a> ah;
    const hash<b> bh;
public:
    hash() : ah(), bh() {}
    size_t operator()(const std::pair<a, b> &p) const {
        std::size_t h=0;
        hash_combine(h,p.first,p.second);
        return h;
    }
};

} // namespaces

class SPHEngine{
public:
    typedef  std::vector<particle> CONTAINER;
    typedef  std::vector<particle*> CONTAINER_P;

    typedef signed long l;
    typedef signed long long ll;

    typedef unsigned long ul;
    typedef unsigned long long ull;
    typedef std::pair<l,std::pair<l,l>> triple;

    std::unordered_map<triple, CONTAINER_P> hashmap;
    std::unordered_map<triple, CONTAINER_P> neighbors_map;

    CONTAINER_P particles;
    std::vector<surface> solid_surfaces;
    double inc=-0.1;
    VEC G;
    double TIME=0.0d;


    double radius; //intorno di influenza (smoothing kernel distance)
    double radius_2;


    double sizeHB; //size of the logical block of the HM

    SPHEngine() : TIME(0.0) {
        G.x=0;
        G.y=-GRAVITY_ACCELERATION*3;
        G.z=0;

        radius=0.0;
        radius_2 = radius*radius;
        sizeHB = 2*radius;



    }

    void setRadius(const double _radius){
        radius=_radius;
        radius_2 = radius*radius;
        sizeHB = 2*radius;
    }







    triple phash(const particle& p){
        triple t;
        double inc = p.pos.x >= 0 ? -0.1 : +0.1;
        double hx = (1+p.pos.x)/(sizeHB);
        double hy = (1+p.pos.y)/(sizeHB);
        double hz = (1+p.pos.z)/(sizeHB);
        t.first = hx;
        t.second.first = hy;
        t.second.second = hz;
        using std::cout; using std::endl;
        // printf("(%f,%f,%f) = >%d:%d:%d\n",p.pos.x,p.pos.y,p.pos.z,t.first, t.second.first,t.second.second);
        return t;
    }

    template<class CONTAINER>
    CONTAINER_P& getNeighbors(particle const& p){

        return particles;
    }

    //naive implementation
    template<class CONTAINER>
    CONTAINER_P& getNeighbors1(particle const& p){
        triple t = phash(p);
        triple t1 = t;
        if(neighbors_map.find(t)==neighbors_map.end()){
            //biuld container
            CONTAINER_P c;
            for (int i = -1; i <= 1; i++) {
                t1.first = t.first+i;
                for (int j = -1; j <= 1; ++j) {
                    t1.second.first = t.second.first +j;
                    for (int k = -1; k <= 1; ++k) {
                        t1.second.second=t.second.second+k;
                        if(hashmap.find(t1) != hashmap.end())
                            c.insert(c.end(),hashmap[t1].begin(),hashmap[t1].end());
                    }
                }
            }
            neighbors_map[t] = c;
        }

        return neighbors_map[t];


        //return hashmap[phash(p)];
        //return particles;
    }

    bool isNeigh(const particle& p, const particle& p1){
        if(glm::distance(p.pos,p1.pos) <= radius)
            return true;
        return false;
    }



    void computeDensity(){

        //  G=glm::rotate(G,1.0,VEC(1.0f, 0.0f, 0.0f));
       // typedef std::vector<particle> CONTAINER;
        for ( auto& p1 : particles) { //O(n)
            auto p = *p1;
            CONTAINER_P neighbors = getNeighbors<CONTAINER_P>(p);  //O(n)

            p.density = 0.0d;
            for ( particle* p1 : neighbors) {
                particle n = *p1;
                if(p.id != n.id && isNeigh(p,n)){
                    VEC diff = p.pos - n.pos;
                    double dist_2 = glm::dot(diff,diff);
                    p.density +=Wpoly6(dist_2, radius);
                }
            }//neigh

            p.density *=p.mass;
            p.pressure =  STIFFNESS * (p.density -  REST_DENSITY);
            //std::cout<<p.density<<std::endl;


        }//particles
    }


    VEC computeExternalFoces(const particle &p){
        VEC _a_extern(0,0,0);
        // std::cout<<std::endl;
        for(auto& w : solid_surfaces){
            double d = glm::dot((w.pos - p.pos),w.norm) +p.pradius; // particle radius

            // std::cout<<w.pos.y<<" "<<p.pos.x<<"= "<<d<<std::endl;
            if(d > 0.0){
                _a_extern +=  WALL_K * w.norm * d;
                _a_extern +=  WALL_DAMPING * glm::dot(p.vel,w.norm) * w.norm;
            }
        }

        return _a_extern;
    }

    void computePressureAcceleration(){
        typedef std::vector<particle> CONTAINER;
        for ( auto& p1 : particles) {
            auto p = *p1;
            //forces
            VEC f_gravity = G * p.density;
            VEC f_pressure=VEC(0,0,0);
            VEC f_viscosity=VEC(0,0,0);
            VEC f_surface=VEC(0,0,0);
            VEC a_external=VEC(0,0,0);
            VEC colorFieldNormal=VEC(0,0,0);
            double colorFieldLaplacian;


            CONTAINER_P neighbors = getNeighbors<CONTAINER_P>(p);

            int cn=0;
            for ( auto& p1 : neighbors) {
                auto n = *p1;
                if(isNeigh(p,n)){
                    cn++;
                    VEC diff = p.pos - n.pos;
                    double dist_2 = glm::dot(diff,diff);

                    VEC poly6Gradient;
                    Wpoly6Gradient(diff, dist_2, poly6Gradient,radius);

                    VEC spikyGradient;
                    WspikyGradient(diff, dist_2, spikyGradient, radius);

                    if(p.id != n.id){
                        f_pressure +=(  p.pressure/pow(p.density,2)+
                                        n.pressure/pow(n.density,2)
                                        )*spikyGradient;

                        f_viscosity += ( n.vel -
                                         p.vel
                                         )* WviscosityLaplacian(dist_2,radius) / n.density;
                    }

                    colorFieldNormal += poly6Gradient / n.density;
                    colorFieldLaplacian += Wpoly6Laplacian(dist_2,radius) / n.density;

                }
            }//neigh


            f_pressure *= -p.mass* p.density;
            f_viscosity *=  VISCOSITY * p.mass;

            colorFieldNormal *= p.mass;
            p.normal = -1.0 * colorFieldNormal;
            colorFieldLaplacian *=p.mass;

            // surface tension force
            double colorFieldNormalMagnitude = colorFieldNormal.length(); //check that this is actually the magnitude of the vector

            if (colorFieldNormalMagnitude >  SURFACE_THRESHOLD) {

                p.flag = true;
                f_surface = - SURFACE_TENSION * colorFieldNormal / colorFieldNormalMagnitude * colorFieldLaplacian;

            }
            else {
                p.flag = false;
            }

            //if(cn <=1)
            //  std::cout<<"NON BUONO"<<std::endl;
            // if(f_pressure.x > 100000){
            //    std::cout<<f_pressure.x<<" "<<f_pressure.y<<" "<<f_pressure.z<<" "<<f_gravity.y<<" "<<p.density<<std::endl;
            // }

            //add up all the forces and compute the resulting acceleration vector
            const double eps = 0.1;
            // if(cn >1 && (p.density >= eps || p.density <= -eps))
            p.acc = (f_pressure  +f_viscosity + f_surface + f_gravity ) / p.density;

             std::cout<<p.acc.x<<" "<<p.acc.y<<std::endl;

            /**add external forces, collision, user interaction, moving rigid bodys etc.**/
            a_external = computeExternalFoces(p);
            p.acc+= a_external;




        }
    }

    void advance(const double dt){
        hashmap.clear();
        for ( auto& p1 : particles) {
            auto p = *p1;
            VEC newPosition = p.pos + p.vel*dt + p.acc*dt*dt;
            VEC newVelocity = (newPosition - p.pos) /dt;

            p.pos = newPosition;
            p.vel = newVelocity;
            if(hashmap.find(phash(p))!= hashmap.end()){
                hashmap[phash(p)].push_back(&p);
            }else{
                CONTAINER_P c(1); c[0] = &p;
                hashmap[phash(p)]=c;
            }
        }
    }

    void step(const double dt){

        computeDensity();
        computePressureAcceleration();
        advance(dt);
        neighbors_map.clear();

    }

    void evolveSTUPID(const int DIMY){
        srand(time(0));
        for ( auto& p1 : particles) {
            auto p = *p1;
            double r = (double)rand()/RAND_MAX*inc;
            if(p.pos.y > DIMY/2)
                inc=-0.1+r;
            if(p.pos.y < -DIMY/2)
                inc=+0.1+r;

            p.pos.y+=inc;
        }
    }

    void addBoundaries(const double DIMX, const double  DIMY, const double  DIMZ){
        double div=4;

        solid_surfaces.push_back(surface(VEC(0,0,1), VEC(0,0,-DIMZ/2))); // back
        solid_surfaces.push_back(surface(VEC(0,0,-1), VEC(0,0,DIMZ/2))); //front

        solid_surfaces.push_back(surface(VEC(1,0,0), VEC(-DIMX/2,0,0))); //left
        solid_surfaces.push_back(surface(VEC(-1,0,0), VEC(DIMX/2,0,0))); //right


        solid_surfaces.push_back(surface(VEC(0,1,0), VEC(0,-DIMY/2,0))); // bottom
        solid_surfaces.push_back(surface(VEC(0,-1,0), VEC(0,DIMY/2,0))); // top

        //  solid_surfaces.push_back(surface(VEC(0,1,0), VEC(0,-DIMY/4,0))); // top

    }

    double frand(double fMin, double fMax)
    {
        double f = (double)rand() / RAND_MAX;
        return fMin + f * (fMax - fMin);
    }
    void initialize(const double DIMX, const double  DIMY, const double  DIMZ, const double MASS,const double PRADIUS){
        srand(time(0));
        addBoundaries(DIMX,DIMY,DIMZ);
        double m,M;

        setRadius(PRADIUS);
        long long id=0;
        double div=2;
        double p = 0.01;



#define ROWS (0.2d)
#define COLS (0.05d )
#define LAYERS (0.021d)
        particle prec;
        for (double y = -COLS; y < COLS; y+= radius/2) {
            for (double x = -ROWS; x < ROWS; x += radius/2) {
                for (double z = -LAYERS; z < LAYERS; z+= radius/2) {
                    m=-PRADIUS/2;
                    M=+PRADIUS/2;
                    particle p;
                    // VEC d(x+frand(m,M), y+frand(m,M),z+5*frand(m,M));
                    VEC d(x,y,z);
                    p.pos=d;
                    p.mass=MASS;
                    p.id=id++;

                    particles.push_back(new particle(p));
                    hashmap[phash(p)].push_back(particles[particles.size()-1]);

                    bool n = isNeigh(prec,p);
                    prec=p;


                }
            }
        }
        using std::cout;
        using std::endl;
        /*  ll c=0;
        for(auto p : hashmap){
            std::cout<<p.second.size()<<std::endl;
            c+=p.second.size();

        }
        std::cout<<c<<endl;*/

        cout << "Loaded cube scenario" << endl;
        cout << "Simulating " << particles.size() << " particles" << endl;

        /* srand(time(0));

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {


                particle p;
                p.pos.x=-DIMX/2+i*f+MASS;
                p.pos.y=DIMY/2+MASS;
                p.pos.z=-DIMZ/2+MASS+j*f;
                p.mass=MASS;
                particles.push_back(p);
            }
        }*/


    }



};

#endif // _SPH_ENGINE
