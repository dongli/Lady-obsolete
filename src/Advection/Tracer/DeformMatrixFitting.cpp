#include "DeformMatrixFitting.h"
#include "TracerSkeleton.h"

namespace lady {
    
#ifdef DEBUG
    double DeformMatrixFitting::maxObjectiveValue = -999;
#endif

DeformMatrixFitting::DeformMatrixFitting(const LADY_DOMAIN &domain,
                                         Tracer *tracer) {
    data.domain = &domain;
    data.x.resize(tracer->getSkeleton().getYs().size());
    data.y.resize(data.x.size());
    for (int i = 0; i < data.x.size(); ++i) {
        data.x[i].set_size(domain.getNumDim());
        data.y[i].set_size(domain.getNumDim());
    }
    // y is fixed
    for (int i = 0; i < data.y.size(); ++i) {
        data.y[i] = (*tracer->getSkeleton().getYs()[i])();
    }
    // create optimization object
    a = new nlopt::opt(nlopt::LD_SLSQP, pow(domain.getNumDim(), 2));
    a->set_min_objective(DeformMatrixFitting::objective, &data);
    a->add_inequality_constraint(DeformMatrixFitting::constraint, &data);
    a->set_stopval(1.0e-6);
    a->set_ftol_rel(1.0e-3);
    a->set_ftol_abs(1.0e-6);
}
    
DeformMatrixFitting::~DeformMatrixFitting() {
    delete a;
}

void DeformMatrixFitting::fit(const TimeLevelIndex<2> &timeIdx,
                              Tracer *tracer) {
    data.count = 0; // iteration counter
    // prepare the first guess
    vector<double> H0(pow(data.domain->getNumDim(), 2));
    int k = 0;
    for (int j = 0; j < data.domain->getNumDim(); ++j) {
        for (int i = 0; i < data.domain->getNumDim(); ++i) {
            H0[k++] = tracer->getH(timeIdx-1)(i, j);
        }
    }
    // prepare data x
    const LADY_SPACE_COORD &x0 = tracer->getX(timeIdx);
    TracerSkeleton &s = tracer->getSkeleton();
    if (dynamic_cast<const geomtk::SphereDomain*>(data.domain) != NULL) {
        for (int i = 0; i < data.x.size(); ++i) {
            data.domain->project(geomtk::SphereDomain::STEREOGRAPHIC, x0,
                                 *s.getXs(timeIdx)[i], data.x[i]);
        }
    } else {
        REPORT_ERROR("Under construction!");
    }
//    if (tracer->getID() == 10000) {
//        data.debug = true;
//        tracer->getH(timeIdx-1).print();
//        std::ofstream file;
//        file.open("x.dat");
//        for (int i = 0; i < data.x.size(); ++i) {
//            file << data.x[i] << endl;
//        }
//    } else {
//        data.debug = false;
//    }
    // perform optimization
    double obj;
    nlopt::result res;
    try {
        res = a->optimize(H0, obj);
    } catch (const std::exception &e) {
        REPORT_ERROR("Encounter exception \"" << e.what() << "\" from NLopt!");
    }
    if (!(res == nlopt::SUCCESS || res == nlopt::FTOL_REACHED ||
          res == nlopt::STOPVAL_REACHED)) {
        REPORT_ERROR("Optimization failed!");
    }
#ifdef DEBUG
    if (maxObjectiveValue < obj) maxObjectiveValue = obj;
#endif
//    if (data.debug) {
//        CHECK_POINT
//    }
    // copy result
    k = 0;
    for (int j = 0; j < data.domain->getNumDim(); ++j) {
        for (int i = 0; i < data.domain->getNumDim(); ++i) {
            tracer->getH(timeIdx)(i, j) = H0[k++];
        }
    }
}

double DeformMatrixFitting::objective(unsigned int n,  const double *x,
                                      double *grad, void *data_) {
    Data *data = static_cast<Data*>(data_);
    mat H(x, data->domain->getNumDim(), data->domain->getNumDim());
    double f = 0.0;
    if (dynamic_cast<const geomtk::SphereDomain*>(data->domain) != NULL) {
        if (data->domain->getNumDim() == 2) {
            // calculate objective value
            for (int i = 0; i < data->x.size(); ++i) {
                    f += pow(norm(data->x[i]-H*data->y[i], 2), 2);
            }
            // calculate gradient
            if (grad != NULL) {
                int k = 0;
                for (int j = 0; j < data->domain->getNumDim(); ++j) {
                    for (int i = 0; i < data->domain->getNumDim(); ++i) {
                        double s0 = 0.0, s1 = 0.0, s2 = 0.0;
                        for (int l = 0; l < data->y.size(); ++l) {
                            s0 += data->y[l](j)*data->y[l](0);
                            s1 += data->y[l](j)*data->y[l](1);
                            s2 += data->y[l](j)*data->x[l](i);
                        }
                        grad[k++] = 2*(s0*H(i, 0)+s1*H(i, 1)-s2);
                    }
                }
            }
        } else if (data->domain->getNumDim() == 3) {
            REPORT_ERROR("Under construction!");
        }
    }
//    if (data->debug) {
//        cout << "[Notice]: Iteration " << ++data->count;
//        cout << ": Objective function value is " << f << endl;
//        // output H into a file for debugging
//        std::ofstream file;
//        if (data->count == 1) {
//            file.open("H.dat");
//        } else {
//            file.open("H.dat", std::ofstream::out | std::ofstream::app);
//        }
//        file << H << endl;
//    }
    return f;
}

double DeformMatrixFitting::constraint(unsigned int n,  const double *x,
                                       double *grad, void *data_) {
    Data *data = static_cast<Data*>(data_);
    mat H(x, data->domain->getNumDim(), data->domain->getNumDim());
    double c = -det(H);
    if (data->domain->getNumDim() == 2) {
        if (grad != NULL) {
            grad[0] =  H(1, 1);
            grad[1] = -H(0, 1);
            grad[2] = -H(1, 0);
            grad[3] =  H(0, 0);
        }
    } else if (data->domain->getNumDim() == 3) {
        if (grad != NULL) {
            grad[0] = H(1, 2)*H(2, 1)-H(1, 1)*H(2, 2);
            grad[1] = H(0, 1)*H(2, 2)-H(0, 2)*H(2, 1);
            grad[2] = H(0, 2)*H(1, 1)-H(0, 1)*H(1, 2);
            grad[3] = H(1, 0)*H(2, 2)-H(1, 2)*H(2, 0);
            grad[4] = H(0, 2)*H(2, 0)-H(0, 0)*H(2, 2);
            grad[5] = H(0, 0)*H(1, 2)-H(0, 2)*H(1, 0);
            grad[6] = H(1, 1)*H(2, 0)-H(1, 0)*H(2, 1);
            grad[7] = H(0, 0)*H(2, 1)-H(0, 1)*H(2, 0);
            grad[8] = H(0, 1)*H(1, 0)-H(0, 0)*H(1, 1);
        }
    }
    return c;
}

}
