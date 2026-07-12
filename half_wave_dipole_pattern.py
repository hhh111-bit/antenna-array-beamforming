"""
Visualize the far-field radiation pattern of a thin half-wave dipole.

Geometry used here:
    - Cartesian coordinates: x, y, z.
    - The dipole is a thin straight wire placed on the z-axis.
    - Its center is at the origin.
    - Its two ends are at z = -lambda/4 and z = +lambda/4.
    - theta is the angle between the observation direction r_hat and
      the positive z-axis, namely the dipole axis.
    - phi is the azimuth angle in the x-y plane.

For a center-fed thin linear dipole of length L, the far-zone field shape is

    F(theta) = [cos((kL/2) cos(theta)) - cos(kL/2)] / sin(theta)

where k = 2*pi/lambda.  For L = lambda/2, kL/2 = pi/2 and cos(kL/2)=0, so

    F(theta) = cos((pi/2) cos(theta)) / sin(theta)

This script plots the normalized magnitude |F(theta)|.  The singular-looking
points at theta = 0 and pi are removable limits; physically the radiation is
zero on the dipole axis.
"""

import numpy as np
import matplotlib.pyplot as plt


def half_wave_dipole_field(theta):
    """Normalized field magnitude of a half-wave dipole."""
    sin_theta = np.sin(theta)
    numerator = np.cos(0.5 * np.pi * np.cos(theta))

    field = np.zeros_like(theta, dtype=float)
    valid = np.abs(sin_theta) > 1e-12
    field[valid] = numerator[valid] / sin_theta[valid]

    return np.abs(field)


def main():
    theta = np.linspace(0.0, np.pi, 361)
    phi = np.linspace(0.0, 2.0 * np.pi, 721)
    theta_grid, phi_grid = np.meshgrid(theta, phi, indexing="ij")

    radius = half_wave_dipole_field(theta_grid)
    radius /= radius.max()

    x = radius * np.sin(theta_grid) * np.cos(phi_grid)
    y = radius * np.sin(theta_grid) * np.sin(phi_grid)
    z = radius * np.cos(theta_grid)

    fig = plt.figure(figsize=(14, 8))

    ax3d = fig.add_subplot(1, 3, 1, projection="3d")
    surface = ax3d.plot_surface(
        x,
        y,
        z,
        facecolors=plt.cm.viridis(radius),
        rstride=3,
        cstride=6,
        linewidth=0,
        antialiased=True,
        shade=False,
    )
    surface.set_facecolor((0, 0, 0, 0))

    # Draw the dipole wire on the z-axis.
    ax3d.plot([0, 0], [0, 0], [-0.45, 0.45], color="black", linewidth=3)
    ax3d.text(0, 0, 0.52, "dipole axis: z", ha="center")

    ax3d.set_title("3D pattern: |F(theta)|")
    ax3d.set_xlabel("x")
    ax3d.set_ylabel("y")
    ax3d.set_zlabel("z")
    ax3d.set_xlim(-1, 1)
    ax3d.set_ylim(-1, 1)
    ax3d.set_zlim(-1, 1)
    if hasattr(ax3d, "set_box_aspect"):
        ax3d.set_box_aspect((1, 1, 1))
    ax3d.view_init(elev=25, azim=35)

    ax_e = fig.add_subplot(1, 3, 2, projection="polar")
    e_cut = half_wave_dipole_field(theta)
    e_cut /= e_cut.max()
    ax_e.plot(theta, e_cut, linewidth=2, label="phi = 0")
    ax_e.plot(theta + np.pi, e_cut, linewidth=2)
    ax_e.set_title("E-plane cut")
    ax_e.set_theta_zero_location("N")
    ax_e.set_theta_direction(-1)
    ax_e.set_rlim(0, 1)
    ax_e.grid(True)

    ax_h = fig.add_subplot(1, 3, 3, projection="polar")
    h_phi = np.linspace(0.0, 2.0 * np.pi, 721)
    h_cut = np.ones_like(h_phi)
    ax_h.plot(h_phi, h_cut, linewidth=2)
    ax_h.set_title("H-plane cut: theta = 90 deg")
    ax_h.set_rlim(0, 1.05)
    ax_h.grid(True)

    fig.suptitle(
        "Thin half-wave dipole, placed along z-axis, centered at origin",
        fontsize=14,
    )
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
