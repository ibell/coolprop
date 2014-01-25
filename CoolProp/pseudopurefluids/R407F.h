
#ifndef R407F_H
#define R407F_H

    class R407FClass : public Fluid{

    public:
        R407FClass();
        ~R407FClass(){};
        double psatL(double);
        double psatV(double);
        double rhosatL(double);
        double rhosatV(double);
    };
#endif
