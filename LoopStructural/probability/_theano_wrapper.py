import theano.tensor as tt
import numpy as np

from ._gradient_calculator import gradients
class LogLikelihood(tt.Op):
    itypes = [tt.dvector]  # expects a vector of parameter values when called
    otypes = [tt.dscalar]  # outputs a single scalar value (normal_likelihood)

    def __init__(self, normal_loglikelihood):
        """
        Initialise the class with the

        Parameters
        ----------
        normal_loglikelihood:
            The log-likelihood function
        model:
            the Geological Model
        sigma:
            The standard deviation
        """
        # add inputs as class attributes
        self.loglike = normal_loglikelihood
        # initialise the LogLikelihoodGradient Op class that calculates the gradient
        # of the our log-likelihood function
        self.loglikegrad = LogLikelihoodGradient(self.loglike)
        self.reps = 1e-3
        self.reltol = 1e-3
        self.epsscale = 0.5
        self.mineps = 1e-9
    def perform(self, node, inputs, outputs):
        # the method that is used when calling the Op
        # (theta,) = inputs
        (theta,) = inputs

        loglikelihood = self.loglike(theta)

        outputs[0][0] = np.array(loglikelihood)  # output the misfit

    def grad(self, inputs, g):
        # the method that calculates the gradients - it actually returns the
        # vector-Jacobian product - g[0] is a vector of parameter values
        (theta,) = inputs  # our parameters
        return [g[0] * self.loglikegrad(theta)]

class LogLikelihoodGradient(tt.Op):
    """
    This Op will be called with a vector of values and also return a vector of
    values - the gradients in each dimension.
    """

    itypes = [tt.dvector]
    otypes = [tt.dvector]

    def __init__(self, normal_loglikelihood):
        """
        Initialise the class with the

        Parameters
        ----------
        normal_loglikelihood:
            The log-likelihood function
        model:
            the Geological Model
        sigma:
            The standard deviation
        """

        # add inputs as class attributes
        self.loglike = normal_loglikelihood


    def perform(self, node, inputs, outputs):
        (theta,) = inputs

        # define version of normal_likelihood function to pass to derivative function
        # values are theta parameters of the normal_likelihood function
        def lnlike(values):
            return self.loglike(values)

        # calculate gradients
        self.reps = 1e-3
        self.reltol = 1e-3
        self.epsscale = 0.5
        self.mineps = 1e-9
        grads = gradients(theta, lnlike)#,reps=self.reps,mineps=self.mineps,epsscale=self.epsscale,reltol=self.reltol)

        outputs[0][0] = grads