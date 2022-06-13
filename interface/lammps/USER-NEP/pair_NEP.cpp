#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <vector>
#include "pair_NEP.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "domain.h"
#include "nep.h"

using namespace LAMMPS_NS;

PairNEP::PairNEP(LAMMPS *lmp) : Pair(lmp)
{
    restartinfo = 0;
    manybody_flag = 1;

    single_enable = 0;

    inited = false;
    allocated = 0;
}

PairNEP::~PairNEP()
{
    if (copymode) return;

    if (allocated) {
        memory->destroy(setflag);
        memory->destroy(cutsq);
    }

}

void PairNEP::allocate()
{
    int n = atom->ntypes;

    memory->create(setflag, n+1, n+1, "pair:setflag");
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= n; j++)
            setflag[i][j] = 1;

    memory->create(cutsq, n+1, n+1, "pair:cutsq");

    allocated = 1;
}

void PairNEP::coeff(int narg, char **arg)
{
    if (!allocated) allocate();
}

void PairNEP::settings(int narg, char **arg)
{
    if (narg != 1) error->all(FLERR, "Illegal pair_style command");
    strcpy(model_filename, arg[0]);
}

void PairNEP::init_style()
{
    int irequest = neighbor->request(this,instance_me);
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->full = 1;
    
    // if (inited) ~nep_model();
    nep_model = NEP3(model_filename);
    inited = true;
    cutoff = nep_model.paramb.rc_radial;
    cutoffsq = cutoff * cutoff;
    int n = atom->ntypes;
    for (int i=1; i<=n; i++)
        for (int j=1; j<=n; j++)
            cutsq[i][j] = cutoffsq;

    if (comm->nprocs != 1) error->all(FLERR, "no parallel plz");
}

double PairNEP::init_one(int i, int j)
{
    return cutoff;
}

void PairNEP::compute(int eflag, int vflag)
{
    if (eflag || vflag) ev_setup(eflag,vflag);
    else evflag = vflag_fdotr = eflag_global = eflag_atom = 0;

    double **x  = atom->x;
    double **f  = atom->f;
    int *atom_type_   = atom->type; 
    int n_atoms = list->inum;

    std::vector<double> lattice(9);
    lattice[0] = domain->xprd; lattice[3] = 0.0;          lattice[6] = 0.0;
    lattice[1] = domain->xy;   lattice[4] = domain->yprd; lattice[7] = 0.0;
    lattice[2] = domain->xz;   lattice[5] = domain->yz;   lattice[8] = domain->zprd;

    std::vector<double> position(n_atoms * 3);
    std::vector<int> atom_type(n_atoms);
    for (int i=0; i<n_atoms; i++)
    {
        position[i]           = x[i][0];
        position[i+n_atoms]   = x[i][1];
        position[i+n_atoms*2] = x[i][2];
        atom_type[i] = atom_type_[i] - 1;
    }

    std::vector<double> nep_energies(n_atoms);
    std::vector<double> nep_forces(n_atoms * 3);
    std::vector<double> nep_virials(n_atoms * 9);

    nep_model.compute(atom_type, lattice, position, nep_energies, nep_forces, nep_virials);

    for (int i=0; i<n_atoms; i++)
    {
        f[i][0] = nep_forces[i];
        f[i][1] = nep_forces[i+n_atoms];
        f[i][2] = nep_forces[i+n_atoms*2];
    }

    if (eflag)
    {
        if (eflag_global)
        {
            double en = 0.0;
            for (int i = 0; i < n_atoms; i++) en += nep_energies[i];
            eng_vdwl = en;
        }

        if (eflag_atom)
        {
            for (int i = 0; i < n_atoms; i++) eatom[i] = nep_energies[i];
        }
    }

    if (vflag)
    {
        if (vflag_global || vflag_fdotr)
        {
            for (int i = 0; i < n_atoms; i++)
            {
                virial[0] += nep_virials[i];
                virial[1] += nep_virials[i+n_atoms*4];
                virial[2] += nep_virials[i+n_atoms*8];
                virial[3] += (nep_virials[i+n_atoms*1] + nep_virials[i+n_atoms*3]) / 2;
                virial[4] += (nep_virials[i+n_atoms*2] + nep_virials[i+n_atoms*6]) / 2;
                virial[5] += (nep_virials[i+n_atoms*5] + nep_virials[i+n_atoms*7]) / 2;
            }
        }

        if (vflag_atom)
        {
            for (int i = 0; i < n_atoms; i++)
            {
                vatom[i][0] = nep_virials[i+n_atoms*0];
                vatom[i][1] = nep_virials[i+n_atoms*4];
                vatom[i][2] = nep_virials[i+n_atoms*8];
                vatom[i][3] = (nep_virials[i+n_atoms*1] + nep_virials[i+n_atoms*3]) / 2;
                vatom[i][4] = (nep_virials[i+n_atoms*2] + nep_virials[i+n_atoms*6]) / 2;
                vatom[i][5] = (nep_virials[i+n_atoms*5] + nep_virials[i+n_atoms*7]) / 2;
            }
        }
    }
}
