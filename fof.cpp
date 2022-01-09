//
// fof.cpp
//
// simple program to implement FoF clustet finder
// with constant link length
//
// Jun Makino Jan 9 2022
//


#include<iostream>
#include<fstream>
#include<unistd.h>
#include<sys/stat.h>
#include<particle_simulator.hpp>


class FileHeader{
public:
    PS::S64 n_body;
    PS::F64 time;
    PS::S32 readAscii(FILE * fp){
        fscanf(fp, "%lf\n", &time);
        fscanf(fp, "%lld\n", &n_body);
        return n_body;
    }
    void writeAscii(FILE* fp) const{
        fprintf(fp, "%e\n", time);
        fprintf(fp, "%lld\n", n_body);
    }
};


class FPFOF{
public:
    PS::S64 id;
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64 search_radius;
    PS::S64 cluster_id;
    PS::S64 cluster_id_prev;
    PS::S64 nneighbours;
    void clear(){
	cluster_id = -1;
	nneighbours=0;
    }
    PS::F64 getRSearch() const{
        return this->search_radius;
    }
    PS::F64vec getPos() const { return pos; }
    void setPos(const PS::F64vec & p) { pos = p; }
    void copyFromForce(const FPFOF & force){
        cluster_id = force.cluster_id;
	nneighbours = force.nneighbours;
    }
    void writeAscii(FILE* fp) const{
	fprintf(fp, "%lld\t%.15g\t%.15g\t%.15g\t%.15g\t%lld\t%lld\n", 
                this->id, this->mass, this->pos.x, this->pos.y, this->pos.z,
		this->nneighbours, this->cluster_id);
    }

    void readAscii(FILE* fp){
	char readbuf[2048];
	fgets(readbuf, 2000, fp);
	sscanf(readbuf, "%lld%lf%lf%lf%lf", 
               &this->id, &this->mass, &this->pos.x, &this->pos.y, &this->pos.z);
    }
    void copyFromFP(const FPFOF & fp){ 
        mass = fp.mass;
        pos = fp.pos;
        id = fp.id;
        cluster_id = fp.cluster_id;
        search_radius = fp.search_radius;
	nneighbours = fp.nneighbours;
    }
    
};

void CalcForceFpFp(const FPFOF * ep_i,
                   const PS::S32 n_ip,
                   const FPFOF * ep_j,
                   const PS::S32 n_jp,
                   FPFOF * force){
    for(PS::S32 i=0; i<n_ip; i++){
        PS::F64vec xi = ep_i[i].pos;
        PS::S64 nn=0;
	
        PS::S64 idminneighbours=ep_i[i].cluster_id;
	PS::F64 r0sq = ep_i[i].search_radius*ep_i[i].search_radius;
        for(PS::S32 j=0; j<n_jp; j++){
            PS::F64vec rij = xi - ep_j[j].pos;
            PS::F64 r2 = rij * rij;
            if (r2 < r0sq){
		nn++;
		if (idminneighbours> ep_j[j].cluster_id){
		    idminneighbours= ep_j[j].cluster_id;
		}
	    }
        }
        force[i].nneighbours += nn;
	if (nn > 0){
	    if (force[i].cluster_id<0 ||force[i].cluster_id > idminneighbours){
		force[i].cluster_id =idminneighbours;
	    }
	}else{
	    if (force[i].cluster_id<0) {
		force[i].cluster_id=ep_i[i].cluster_id;
	    }
	}

    }
}

template<class Tpsys>
void initialize_cluster_id(Tpsys & system)
{
    PS::S32 n = system.getNumberOfParticleLocal();
    for(int i=0; i<n; i++){
        system[i].cluster_id = system[i].id;
        system[i].cluster_id_prev =system[i].id;
    }
}

template<class Tpsys>
void set_rsearch(Tpsys & system,
		 PS::F64 rsearch)
{
    PS::S32 n = system.getNumberOfParticleLocal();
    for(int i=0; i<n; i++){
        system[i].search_radius = rsearch;
    }
}

template<class Tpsys>
void copy_cluster_id(Tpsys & system)
{
    PS::S32 n = system.getNumberOfParticleLocal();
    for(int i=0; i<n; i++){
        system[i].cluster_id_prev =  system[i].cluster_id;
    }
}

template<class Tpsys>
int cluster_id_updated(Tpsys & system)
{
    PS::S32 n = system.getNumberOfParticleLocal();
    int nupdated=0;
    for(int i=0; i<n; i++){
        if (system[i].cluster_id_prev !=  system[i].cluster_id){
	    nupdated++;
	}
    }
    int ntotalupdated= PS::Comm::getSum(nupdated);
    return ntotalupdated;
}



int main(int argc, char *argv[]){
    PS::F64 theta = 0.5;
    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);
    PS::Initialize(argc, argv);
    const PS::S32 n_leaf_limit = 8;
    PS::S32 n_group_limit = 64;
    char sinput[1024];
    char soutput[1024];
    PS::S64 n_loc;
    int c;
    PS::F64 rsearch = 0.01;
    while((c=getopt(argc,argv,"i:o:r:h")) != -1){
        switch(c){
        case 'i':
            strncpy(sinput,optarg,1000);
            break;
        case 'o':
            strncpy(soutput,optarg,1000);
            break;
        case 'r':
            rsearch = atof(optarg);
            std::cerr<<"rsearch="<<rsearch<<std::endl;
            break;
         case 'h':
            std::cerr<<"i: input file name"<<std::endl;
            std::cerr<<"o: output file name"<<std::endl;
            std::cerr<<"r: rsearch (default: 0.01)"<<std::endl;
            std::cerr<<"h: Print this help message"<<std::endl;
            return 0;
        }
    }

    PS::ParticleSystem<FPFOF> system_grav;
    system_grav.initialize();
    PS::S32 n_grav_loc;

    FileHeader header;
    system_grav.readParticleAscii(sinput, header);
    n_loc = system_grav.getNumberOfParticleLocal();
    initialize_cluster_id(system_grav);
    set_rsearch(system_grav,rsearch);


    const PS::F64 coef_ema = 0.3;
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);
    dinfo.collectSampleParticle(system_grav);
    dinfo.decomposeDomain();
    system_grav.exchangeParticle(dinfo);
    n_grav_loc = system_grav.getNumberOfParticleLocal();
    PS::TreeForForceShort<FPFOF, FPFOF, FPFOF>::Scatter tree_grav;
    tree_grav.initialize(n_grav_loc, theta, n_leaf_limit, n_group_limit);
    int nupdated;
    int nloop=0;
    do{
	tree_grav.calcForceAllAndWriteBack(CalcForceFpFp,  system_grav, dinfo);
	nupdated = cluster_id_updated(system_grav);
	fprintf(stderr, "loop %d updated=%d\n", nloop, nupdated);
	copy_cluster_id(system_grav);
	nloop++;
    }while(nupdated);

    system_grav.writeParticleAscii(soutput, header);
    return 0;
}
