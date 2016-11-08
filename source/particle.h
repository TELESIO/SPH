#ifndef _SPH_PARTICLE
#define _SPH_PARTICLE

#include<glm/glm.hpp>


class particle{

public:
    typedef glm::tvec3<double> VEC;
    long long id;
    particle(): pos(0,0,0),vel(0.0,0.0,0.0),acc(0.0,0.0,0.0),normal(0.0,0.0,0.0){
    };

    VEC pos;
    VEC vel;
    VEC acc;
    VEC normal;
    double mass=0; //KG
    double pradius=0;
    double density=0;
    double pressure=0;
    bool flag=false;


};


#endif //_SPH_PARTICLE
