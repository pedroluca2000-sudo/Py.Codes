# eos_helmholtz_short_forms.py
from __future__ import annotations
from dataclasses import dataclass
import math
from typing import List
from eos_helmholtz_core import ResidualTerm, Alpha0Ideal

# ----------------------------------------------------------
# Termos residuais da forma "short" (Span & Wagner / Lemmon & Span)
# Há dois grupos usuais: polinomiais em δ^d τ^t e "exponenciais" com
# fatores tipo exp(-α (δ-ε)^2 - β (τ-γ)^2) e multiplicadores δ^d τ^t.
# ----------------------------------------------------------

@dataclass
class PolyTerm(ResidualTerm):
    n: float
    d: float
    t: float
    c: float  # pode ser 0 na forma "short" sem log(δ)

    def val(self, delta, tau):         return self.n * (delta**self.d) * (tau**self.t)
    def d_delta(self, delta, tau):      return self.n * self.d * (delta**(self.d-1)) * (tau**self.t)
    def d_tau(self, delta, tau):        return self.n * (delta**self.d) * self.t * (tau**(self.t-1))
    def d_delta2(self, delta, tau):     return self.n * self.d*(self.d-1) * (delta**(self.d-2)) * (tau**self.t)
    def d_deltatau(self, delta, tau):   return self.n * self.d * (delta**(self.d-1)) * self.t * (tau**(self.t-1))
    def d_tau2(self, delta, tau):       return self.n * (delta**self.d) * self.t*(self.t-1) * (tau**(self.t-2))

@dataclass
class ExpoTerm(ResidualTerm):
    n: float
    d: float
    t: float
    alpha: float
    beta: float
    gamma: float
    epsilon: float

    def _theta(self, delta, tau):
        # “gaussiano” deslocado típico das formas short
        return math.exp(-self.alpha*(delta - self.epsilon)**2 - self.beta*(tau - self.gamma)**2)

    def val(self, delta, tau):
        return self.n * (delta**self.d) * (tau**self.t) * self._theta(delta, tau)

    def d_delta(self, delta, tau):
        theta = self._theta(delta, tau)
        pref  = self.n * (tau**self.t) * (delta**self.d)
        return pref * theta * ( self.d/delta - 2.0*self.alpha*(delta - self.epsilon) )

    def d_tau(self, delta, tau):
        theta = self._theta(delta, tau)
        pref  = self.n * (delta**self.d) * (tau**self.t)
        return pref * theta * ( self.t/tau - 2.0*self.beta*(tau - self.gamma) )

    def d_delta2(self, delta, tau):
        # derivada em δ duas vezes (produto + cadeia)
        theta = self._theta(delta, tau)
        A = self.d/delta - 2.0*self.alpha*(delta - self.epsilon)
        B = - self.d/(delta**2) - 2.0*self.alpha + 4.0*(self.alpha**2)*(delta - self.epsilon)**2
        pref = self.n * (tau**self.t) * (delta**self.d)
        return pref*theta*(A**2 + B)

    def d_deltatau(self, delta, tau):
        theta = self._theta(delta, tau)
        A = self.d/delta - 2.0*self.alpha*(delta - self.epsilon)
        C = self.t/tau - 2.0*self.beta*(tau - self.gamma)
        pref = self.n * (delta**self.d) * (tau**self.t)
        # derivada mista de produto (ver notas no artigo)
        return pref*theta*(A*C - 2.0*self.alpha*(delta - self.epsilon)*C - 2.0*self.beta*(tau - self.gamma)*A)

    def d_tau2(self, delta, tau):
        theta = self._theta(delta, tau)
        C = self.t/tau - 2.0*self.beta*(tau - self.gamma)
        D = - self.t/(tau**2) - 2.0*self.beta + 4.0*(self.beta**2)*(tau - self.gamma)**2
        pref = self.n * (delta**self.d) * (tau**self.t)
        return pref*theta*(C**2 + D)

# ----------------------------------------------------------
# α0(τ) (ideal-gás) – aqui deixamos uma classe-abstrata a ser
# concretizada no arquivo de parâmetros do fluido, pois cada
# correlação pode optar por uma expressão (p.ex., termos
# logarítmicos + somatórios racionais do tipo ∑ a_i ln(1 - e^{-θ_i τ})).
# ----------------------------------------------------------

class Alpha0_ShortPolar(Alpha0Ideal):
    """
    Implementação típica:
      α0(δ,τ) = ln(δ) + a + b*τ + ∑ c_i ln(1 - exp(-θ_i τ))
    (aqui, o conjunto concreto fica no arquivo de parâmetros)
    """
    def __init__(self, has_ln_delta: bool, a0_eval_func):
        self.has_ln_delta = has_ln_delta
        self._eval = a0_eval_func

    def eval_all(self, delta: float, tau: float):
        a0, a0_tau = self._eval(tau)  # só τ
        a0_delta = 1.0/delta if self.has_ln_delta else 0.0
        return a0 + (math.log(delta) if self.has_ln_delta else 0.0), a0_delta, a0_tau