#ifndef _SPH_SOLID_SURFACE
#define _SPH_SOLID_SURFACE

class surface{
public:
    typedef particle::VEC VEC;
    surface(const VEC& nor, const VEC& po){
        pos=po;
        norm = nor;
    }


    VEC pos;
    VEC norm;

};

#endif//_SPH_SOLID_SURFACE

