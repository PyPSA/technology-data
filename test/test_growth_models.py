# SPDX-FileCopyrightText: 2025 The technology-data authors
# SPDX-License-Identifier: MIT
"""Unit tests for growth models in technologydata.technologies.growth_models."""

import numpy as np
import pytest

from technologydata.technologies.growth_models import (
    ExponentialGrowth,
    GeneralLogisticGrowth,
    GompertzGrowth,
    LinearGrowth,
    LogisticGrowth,
)


def test_linear_add_data_points() -> None:
    """Test adding data points to LinearGrowth model."""
    x0 = 0
    m = 2.0
    A = 1.0

    x = np.arange(-10, 30)
    y = m * x + A

    model = LinearGrowth()
    for xi, yi in zip(x, y):
        model.add_data((xi, yi))
    model.fit(p0={"x0": x0, "m": m, "A": A})

    assert model.x0 == pytest.approx(x0, rel=1e-2)
    assert model.m == pytest.approx(m, rel=1e-2)
    assert model.A == pytest.approx(A, rel=1e-2)


def test_linear_growth_projection() -> None:
    """Test the projection method of LinearGrowth."""
    x0 = 0
    m = 2.0
    A = 5.0
    x = 10
    model = LinearGrowth(x0=x0, m=2.0, A=5.0, data_points=[])
    assert model.project(10) == pytest.approx(m * x + A)


def test_linear_growth_fit() -> None:
    """Test fitting of the LinearGrowth model."""
    x0 = 0
    m = 2.0
    A = 1.0

    x = np.arange(-10, 30)
    y = m * x + A

    model = LinearGrowth(x0=x0, m=None, A=None, data_points=[*zip(x, y)])
    model.fit(p0={"x0": x0, "m": m, "A": A})
    assert model.x0 == pytest.approx(x0, rel=1e-2)
    assert model.m == pytest.approx(m, rel=1e-2)
    assert model.A == pytest.approx(A, rel=1e-2)


def test_exponential_growth_projection() -> None:
    """Test the projection method of ExponentialGrowth."""
    A = 2.0
    m = 1.0
    k = 0.5
    x0 = 0.0
    x = 2

    model = ExponentialGrowth(A=A, m=m, k=k, x0=x0, data_points=[])
    assert model.project(2) == pytest.approx(A + m * np.exp(k * (x - x0)))


def test_exponential_growth_fit() -> None:
    """Test fitting of the ExponentialGrowth model."""
    A = 2.0
    m = 1.0
    k = 0.5
    x0 = 0.0

    x = np.arange(0, 30)
    y = A + m * np.exp(k * (x - x0))

    model = ExponentialGrowth(A=None, m=None, k=None, x0=None, data_points=[*zip(x, y)])
    model.fit(p0={"A": A, "m": m, "k": k, "x0": x0})

    assert model.A == pytest.approx(A, rel=1e-2)
    assert model.m == pytest.approx(m, rel=1e-2)
    assert model.k == pytest.approx(k, rel=1e-2)
    assert model.x0 == pytest.approx(x0, rel=1e-2)


def test_logistic_growth_projection() -> None:
    """Test the projection method of LogisticGrowth."""
    A = 2.0
    L = 10.0
    k = 1.0
    x0 = 0.0

    model = LogisticGrowth(A=A, L=L, k=k, x0=x0, data_points=[])
    assert model.project(0) == pytest.approx(A + L / (1 + np.exp(-k * (0 - x0))))


def test_logistic_growth_fit() -> None:
    """Test fitting of the LogisticGrowth model."""
    A = 2.0
    L = 10.0
    k = 1.0
    x0 = 0.0

    x = np.arange(-2, 30)
    y = A + L / (1 + np.exp(-k * (x - x0)))

    model = LogisticGrowth(A=None, L=None, k=None, x0=None, data_points=[*zip(x, y)])
    model.fit(p0={"A": A, "L": L, "k": k, "x0": x0})

    assert model.A == pytest.approx(A, rel=1e-2)
    assert model.L == pytest.approx(L, rel=1e-2)
    assert model.k == pytest.approx(k, rel=1e-2)
    assert model.x0 == pytest.approx(x0, rel=1e-2)


def test_gompertz_growth_projection() -> None:
    """Test the projection method of GompertzGrowth."""
    A = 10.0
    k = 1.0
    x0 = 0.0
    b = 1.0

    model = GompertzGrowth(A=A, k=k, x0=x0, b=b, data_points=[])

    assert model.project(0) == pytest.approx(A * np.exp(-b * np.exp(-k * (0 - x0))))


def test_gompertz_growth_fit() -> None:
    """Test fitting of the GompertzGrowth model."""
    A = 10.0
    k = 1.0
    x0 = 0.0
    b = 1.0

    x = np.arange(-2, 30)
    y = A * np.exp(-b * np.exp(-k * (x - x0)))

    model = GompertzGrowth(A=None, k=None, x0=None, b=None, data_points=[*zip(x, y)])
    model.fit(p0={"A": A, "k": k, "x0": x0, "b": b})

    assert model.A == pytest.approx(A, rel=1e-2)
    assert model.k == pytest.approx(k, rel=1e-2)
    assert model.x0 == pytest.approx(x0, rel=1e-2)
    assert model.b == pytest.approx(b, rel=1e-2)


def test_general_logistic_growth_projection() -> None:
    """Test the projection method of GeneralLogisticGrowth."""
    x0 = 0
    A = 2.0
    K = 10.0
    B = 1.0
    nu = 1.0
    Q = 1.0
    C = 1.0

    model = GeneralLogisticGrowth(x0=x0, A=A, K=K, B=B, nu=nu, Q=Q, C=C, data_points=[])

    assert model.project(0) == pytest.approx(model.function(0, x0, A, K, B, nu, Q, C))


def test_general_logistic_growth_fit() -> None:
    """Test fitting of the GeneralLogisticGrowth model."""
    x0 = 0
    A = 2.0
    K = 10.0
    B = 1.0
    nu = 1.0
    Q = 1.0
    C = 1.0

    x: np.array = np.arange(-2, 30)
    y: np.array = GeneralLogisticGrowth().function(x, x0, A, K, B, Q, C, nu)

    model = GeneralLogisticGrowth(data_points=[*zip(x, y)])
    model.fit(p0={"x0": x0, "A": A, "K": K, "B": B, "nu": nu, "Q": Q, "C": C})

    assert model.A == pytest.approx(A, rel=1e-2)
    assert model.K == pytest.approx(K, rel=1e-2)
    assert model.B == pytest.approx(B, rel=1e-2)
    assert model.nu == pytest.approx(nu, rel=1e-2)
    assert model.Q == pytest.approx(Q, rel=1e-2)
    assert model.C == pytest.approx(C, rel=1e-2)


# TODO: Add tests for
# - partial fitting (some parameters fixed)
# - fitting with noise/without p0
# - illegal arguments to p0


def test_partial_fitting_linear_growth() -> None:
    """Test fitting LinearGrowth with one parameter fixed and one to fit."""
    x0 = 0
    m = 2.0
    A = 1.0
    x: np.array = np.arange(-10, 30)
    y: np.array = LinearGrowth().function(x, x0, m, A)

    # Fix x0 and m, only fit A
    model = LinearGrowth(x0=x0, m=m, A=None, data_points=[*zip(x, y)])
    model.fit(p0={"A": A * 0.5})

    assert model.m == pytest.approx(m, rel=1e-2)
    assert model.A == pytest.approx(A, rel=1e-2)


def test_fit_with_noisy_data() -> None:
    """Test fitting LinearGrowth with noisy data to check robustness."""
    rng = np.random.default_rng(42)
    x0 = 0
    m = 2.0
    A = 20.0
    x: np.array = np.arange(-10, 30)
    noise = rng.normal(0, 0.1, size=x.shape)
    y: np.array = LinearGrowth().function(x, x0, m, A) + noise

    # Fitting with x0 and A free can lead to competing solutions, so fix x0
    model = LinearGrowth(x0=x0, data_points=[*zip(x, y)])
    model.fit()

    # Should be close to true values despite noise
    assert model.m == pytest.approx(m, rel=1e-1)
    assert model.A == pytest.approx(A, rel=1e-1)


def test_fit_without_p0() -> None:
    """Test fitting LinearGrowth without providing p0 (should use defaults and still work)."""
    x0 = 0
    m = 10.0
    A = 5.0
    x: np.array = np.arange(-25, 25)
    y: np.array = LinearGrowth().function(x, x0, m, A)

    # Fix x0 to avoid competing solutions
    model = LinearGrowth(x0=x0, data_points=[*zip(x, y)])
    model.fit()  # No p0 provided

    assert model.m == pytest.approx(m, rel=1e-2)
    assert model.A == pytest.approx(A, rel=1e-2)


def test_illegal_p0_argument() -> None:
    """Test that fitting works with missing p0 keys (should use defaults, not raise)."""
    x0 = 5
    m = 2.0
    A = 1.0
    x: np.array = np.arange(-10, 30)
    y: np.array = LinearGrowth().function(x, x0, m, A)

    model = LinearGrowth(data_points=[*zip(x, y)])
    # Provide p0 missing one parameter (A missing)
    model.fit(p0={"x0": x0, "m": m})

    # Should still fit correctly
    assert model.m == pytest.approx(m, rel=1e-2)
    assert model.A == pytest.approx(A, rel=1e-2)


def test_fit_with_p0_far_from_target() -> None:
    """Test fitting LinearGrowth with p0 far from the true values (should still converge)."""
    x0 = 0
    m = 20.0
    A = 5.0
    x: np.array = np.arange(-10, 30)
    y: np.array = LinearGrowth().function(x, x0, m, A)

    # Fix x0 to avoid competing solutions
    model = LinearGrowth(x0=x0, data_points=[*zip(x, y)])
    # Provide p0 far from the actual values
    model.fit(p0={"m": m * 10, "A": A * -2})

    assert model.m == pytest.approx(m, rel=1e-2)
    assert model.A == pytest.approx(A, rel=1e-2)
