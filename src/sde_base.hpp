#pragma once


struct SDE {
    virtual double drift(double x, double t) const = 0;
    virtual double diffusion(double x, double t) const = 0;

    virtual double diffusion_derivative(double x, double t) const { return 0.0; }

    virtual ~SDE() = default;
};
