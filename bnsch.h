

void initialization(float ***phi, int sphere_id);

void cahn(float ***c_old, float ***c_new, float ***o1, float ***o2, float ***o3, float BB, float ***ct);

void source(float ***c_old, float ***o1, float ***o2, float ***o3, float ***src_c, float ***src_mu);


void vcycle(float ***uf_new, float ***uf_old, float ***wf_new, 
        float ***su, float ***sw, int nxf, int nyf, int nzf, int ilevel, float BB);


void relax(float ***c_new, float ***c_old, float ***mu_new, float ***su, float ***sw, 
        int ilevel, int nxt, int nyt, int nzt, float BB);

void restrict2(float ***uf, float ***uc, float ***vf, float ***vc, 
        int nxc, int nyc, int nzc);

void restrict3(float ***uf, float ***uc, float ***vf, float ***vc, float ***wf, float ***wc,
        int nxc, int nyc, int nzc);
void prolong_ch(float ***uc, float ***uf, float ***vc, float ***vf, 
             int nxc, int nyc, int nzc);


void defect(float ***duc, float ***dwc,
        float ***uf_new, float ***uf_old, float ***wf_new,
        float ***suf, float ***swf, int nxf, int nyf, int nzf,
        float ***uc_new, float ***uc_old, float ***wc_new,
        int nxc, int nyc, int nzc, float BB);


void nonL(float ***ru, float ***rw, float ***c_new, float ***c_old, float ***mu_new, 
        int nxt, int nyt, int nzt, float BB);

void laplace_ch(float ***a, float ***lap_a, int nxt, int nyt, int nzt);

float error(float ***c_old, float ***c_new, int nxt, int nyt, int nzt);

void augmenc(float ***c, int nxt, int nyt, int nzt);

float Bphi(float ***oc);


