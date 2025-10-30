import numpy as np

from python.validation import test_dhcp_solver as dhcp


def test_small_forward_run_shapes_and_finiteness():
    nt = 3
    ni, nj, nk = 4, 3, 2

    dx = 0.12e-3
    dy = 0.12e-3 * np.sin(np.deg2rad(80.0)) / np.sin(np.deg2rad(45.0))
    dz, dz_b, dz_t = dhcp.build_z_grid(nk, 0.5e-3, 3.0)

    T_init = np.full((ni, nj, nk), 573.15, dtype=np.float64)
    q_surface = np.zeros((nt - 1, ni, nj), dtype=np.float64)

    T_result, diagnostics = dhcp.multiple_time_step_solver_dhcp(
        T_init,
        q_surface,
        nt,
        dhcp.RHO,
        dhcp.CP_COEFFS,
        dhcp.K_COEFFS,
        dx,
        dy,
        dz,
        dz_b,
        dz_t,
        dt=1.0e-3,
        rtol=1.0e-5,
        maxiter=200,
        verbose=False,
    )

    assert T_result.shape == (nt, ni, nj, nk)
    assert len(diagnostics.iterations) == nt - 1
    assert np.isfinite(np.asarray(T_result)).all()
