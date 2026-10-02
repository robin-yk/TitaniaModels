"""Build-time LaTeX typesetting. Native MathML keeps the page self-contained."""
import re
from latex2mathml.converter import convert

EQUATIONS = [
 [r'\mathrm{CO_2}+\mathrm{H_2}\rightleftharpoons\mathrm{CO}+\mathrm{H_2O}'],
 [r'G=\sum_{i\in\mathrm{gas}}n_i\mu_i+\sum_{j\in\mathrm{solid}}m_j\mu_j^\circ(T)',r'\mu_i=\mu_i^\circ(T)+RT\ln\left(\frac{y_iP}{P^\circ}\right)'],
 [r'\mu_i=\sum_e a_{ie}\lambda_e\quad(i\in\mathrm{gas})',r'\mu_j^\circ=t_j\lambda_{\mathrm{Ti}}+o_j\lambda_{\mathrm O}\quad(j\in\mathrm{active\ solids})'],
 [r'n_{\mathrm{V,tot}}=n_{\mathrm{iso}}+n_{\mathrm{assoc}}+n_{\mathrm{below}}'],
 [r'\frac{dq}{dt}=g\left(1-\frac{q}{n_s}\right)',r'q(t)=n_s\left[1-\exp\left(-\frac{gt}{n_s}\right)\right]',r'n_{\mathrm{below}}(t)=gt-q(t)'],
 [r'\frac{dn_{\mathrm{iso}}}{dt}=g\exp\left(-\frac{gt}{n_s}\right)-\frac{2k_{\mathrm{loss}}}{n_s}n_{\mathrm{iso}}^2',r'n_{\mathrm{iso}}(0)=0',r'k_{\mathrm{loss}}=\nu_{\mathrm{eff}}\exp\left(-\frac{E_{\mathrm{loss}}}{k_{\mathrm B}T}\right)',r'n_{\mathrm{assoc}}=q-n_{\mathrm{iso}}'],
 [r'\frac{r_{\mathrm{CO}}(\mathrm{sample})}{r_{\mathrm{CO}}(\mathrm{R600})}=\frac{n_{\mathrm{iso}}(\mathrm{sample})}{n_{\mathrm{iso}}(\mathrm{R600})}'],
]

def block(lines):
    return '<div class="eq">'+''.join(convert(s, display='block') for s in lines)+'</div>'

def typeset(html):
    matches = list(re.finditer(r'<div class="eq">.*?</div>', html, re.S))
    if len(matches) != len(EQUATIONS):
        raise ValueError('Update the display-equation registry after changing the template')
    for match, lines in reversed(list(zip(matches, EQUATIONS))):
        html = html[:match.start()] + block(lines) + html[match.end():]
    html = re.sub(r'\\\[(.*?)\\\]', lambda m: block(m[1].split(';;')), html, flags=re.S)
    html = re.sub(r'\\\((.*?)\\\)', lambda m: convert(m[1]), html)
    return html
