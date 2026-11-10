# User-supplied GAP-port screenshots

`gap-ports-overview.png` shows the instrumented small taper fixture, with internal ports 3 through 7 at the returns of T3, T2, T1, the regular resonator and the output-end resonator. External ports 1 and 2 are retained at the ends. This is not the ten-unit-cell full filter. The two output-side resonator structures belong to the extraction fixture; do not caption the image as the final full-array geometry.

`gap-port-detail.png` shows port 6 at the regular resonator's lower return. The screenshot displays a minus mark above and a plus mark below, defining the port's local reference orientation. These signs are voltage/current reference conventions, not positive and negative physical charges, and do not indicate that a ground plane is held at a positive voltage. The native XML identifies GAP port 6 at (3376,1064) micrometres. The analysis applies the recorded diagonal sign transform [1,1,-1,-1,-1,-1,-1] to keep internal coordinates consistent. Do not infer a separate common ground for all delta-gap ports from this image.

Suggested caption: **Temporary internal electrical access for network characterization.** Delta-gap ports are inserted at shared metal edges without a visible finite gap. Their five interfaces are mathematically shorted when calculating the external response. Their benefit over ordinary two-port sensitivity optimization was not established; the final full-filter simulation has only its two external ports.

Images were supplied through mailbox request 0014 and are preserved unchanged. Canonical native instrumented project: `../local_model/project.sonx`.
