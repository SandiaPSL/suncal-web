''' Custom ipywidgets for Suncal web interface '''

import ipywidgets as widgets
import matplotlib.pyplot as plt
from scipy import stats

import suncal


# Note: LaTeX is not being rendered by MathJax in ipywidgets.Output
# when run in Voila: https://github.com/voila-dashboards/voila/issues/833.
# Results would look better with mathfmt = 'latex', but dosen't work in Voila yet.
reportformat = {'mathfmt': 'svg'}


class FunctionWidget(widgets.HBox):
    def __init__(self, fname='f', fexpr='a + b'):
        self.fname = widgets.Text(fname, layout={'max_width': '50px'})
        self.fexpr = widgets.Text(fexpr)
        self.units = widgets.Text('[units]', layout={'max_width': '150px'})
        super().__init__([self.fname, self.fexpr, self.units])


class FunctionList(widgets.VBox):
    def __init__(self):
        self.funclist = [FunctionWidget()]
        self.btnadd = widgets.Button(description='+', layout={'width': 'max-content'})
        self.btnrem = widgets.Button(description='–', layout={'width': 'max-content'})
        self.btnadd.on_click(self.additem)
        self.btnrem.on_click(self.remitem)
        self.btnnext = widgets.Button(description='Next')
        self.btnlayout = widgets.HBox([self.btnadd, self.btnrem])
        super().__init__([self.btnlayout] + self.funclist + [self.btnnext])
    def additem(self, z):
        item = FunctionWidget()
        self.funclist.append(item)
        self.children = [self.btnlayout] + self.funclist + [self.btnnext]
    
    def remitem(self, z):
        if len(self.funclist) > 1:
            self.funclist.pop(-1)
            self.children = [self.btnlayout] + self.funclist + [self.btnnext]


class UncDistributionWidget(widgets.VBox):
    ''' Distribution Edit Widget for Uncertainty Propagation. Unlike the DIstributionWidget
        used for Risk, this one hides the mean value, and defines normal/t using k and/or confidence.
        Also always shows deg.f and units.
    '''
    def __init__(self, name='a'):
        self.name = widgets.Text(f'u({name})', description='Name')
        self.dist = widgets.Dropdown(options=['normal', 't', 'uniform', 'gamma', 'triangular'], value='normal', description='Distribution')
        self.param1 = widgets.FloatText(value=1, description='Uncertainty', min=0)
        self.param2 = widgets.FloatText(value=2, description='k', step=.01, min=1)
        self.param3 = widgets.FloatText(value=95.45, description='Confidence', step=.01, min=.1)
        self.degf = widgets.FloatText(value=999, min=1, description='Deg. F.')
        self.units = widgets.Text(description='Units')
        self.dist.observe(self.on_distchange, names='value')
        self.param2.observe(self.on_confchange, names='value')
        self.param3.observe(self.on_confchange, names='value')
        self.degf.observe(self.on_confchange, names='value')
        super().__init__([self.name, self.dist, self.param1, self.param2, self.param3, self.degf, self.units])

    def on_confchange(self, change):
        if change['owner'].description == 'k':
            self.param3.unobserve(self.on_confchange, names='value')
            self.param3.value = suncal.ttable.confidence(change['new'], self.degf.value) * 100
            self.param3.observe(self.on_confchange, names='value')
        elif change['owner'].description == 'Confidence':
            self.param2.unobserve(self.on_confchange, names='value')
            self.param2.value = suncal.ttable.t_factor(change['new']/100, self.degf.value)
            self.param2.observe(self.on_confchange, names='value')
        elif change['owner'].description == 'Deg. F.' and self.param2.description == 'k':
            self.param3.value = suncal.ttable.confidence(self.param2.value, change['new']) * 100

    def on_distchange(self, change):
        params = {'normal': ('Uncertainty', 'k', 'Confidence'),
                  't': ('Uncertainty', 'k', 'Confidence'),
                  'gamma': ('a', 'b'),
                  'triangular': ('a',),
                  'uniform': ('a',),
                 }.get(self.dist.value)

        self.param1.description = params[0]
        if len(params) > 2:
            self.param2.description = params[1]
            self.param2.layout.visibility = 'visible'
            self.param3.description = params[2]
            self.param3.layout.visibility = 'visible'
        elif len(params) > 1:
            self.param2.description = params[1]
            self.param2.layout.visibility = 'visible'
            self.param3.layout.visibility = 'hidden'
        else:
            self.param2.layout.visibility = 'hidden'
            self.param3.layout.visibility = 'hidden'

    def get_distargs(self):
        distargs = {'dist': self.dist.value,
                    'name': self.name.value,
                    'units': self.units.value,
                    'df': self.degf.value}

        if self.dist.value in ['normal', 't']:
            distargs.update({'unc': self.param1.value,
                             'k': self.param2.value,
                             'conf': self.param3.value/100})
        else:
            distargs.update({self.param1.description: self.param1.value,
                             self.param2.description: self.param2.value})
        return distargs


class InputValue(widgets.VBox):
    ''' Measured value inputs for one variable '''
    def __init__(self, name='a'):
        self.name = name
        self.nominal = widgets.FloatText(value=1, description='Nominal')
        self.units = widgets.Text(description='Units', continuous_update=False)
        self.uncerts = widgets.Tab(children=[UncDistributionWidget(name)])
        self.uncerts.set_title(0, 'Uncertainty 1')

        self.btnadd = widgets.Button(description='+', layout={'width': 'max-content'})
        self.btnrem = widgets.Button(description='–', layout={'width': 'max-content'})
        self.btnadd.on_click(self.additem)
        self.btnrem.on_click(self.remitem)
        self.btnlayout = widgets.HBox([self.btnadd, self.btnrem])

        self.units.observe(self.on_unitchange, names='value')
        super().__init__([self.nominal, self.units, self.btnlayout, self.uncerts])

    def on_unitchange(self, change):
        units = self.units.value
        for child in self.uncerts.children:
            child.units.value = units

    def additem(self, change):
        n = len(self.uncerts.children)
        self.uncerts.children = list(self.uncerts.children) + [UncDistributionWidget(name=str(n+1))]
        self.uncerts.set_title(n, 'Uncertainty '+str(n+1))
        self.uncerts.selected_index = n

    def remitem(self, change):
        n = self.uncerts.selected_index
        children = list(self.uncerts.children)
        children.pop(n)
        self.uncerts.children = children


class InputList(widgets.VBox):
    def __init__(self):
        self.inputs = []
        self.tab = widgets.Tab()
        self.btnnext = widgets.Button(description='Next')
        self.set_names('a')
        super().__init__([self.tab, self.btnnext])

    def get_names(self):
        return [self.tab.get_title(i) for i in range(len(self.tab.children))]

    def set_names(self, names):
        self.inputs = [InputValue(n) for n in names]
        self.tab.children = self.inputs
        [self.tab.set_title(i, n) for i, n in enumerate(names)]

    def add_inpt(self, name):
        self.inputs.append(InputValue(name))
        self.tab.children = self.inputs
        self.tab.set_title(len(self.tab.children)-1, name)

    def rem_inpt(self, name):
        i = self.get_names().index(name)
        self.inputs.pop(i)
        self.tab.children = self.inputs


class Correlation(widgets.HBox):
    def __init__(self, names):
        self.v1 = widgets.Dropdown(options=names, layout={'max_width': '70px'})
        self.v2 = widgets.Dropdown(options=names, layout={'max_width': '70px'})
        if len(names) > 1:
            self.v2.value = names[1]
        self.val = widgets.FloatText(min=-1, max=1, step=.001, layout={'max_width': '70px'})
        super().__init__([self.v1, self.v2, self.val])


class Settings(widgets.VBox):
    def __init__(self):
        self.names = []
        self.nsamples = widgets.IntText(min=1, max=10000000, value=1000000, description='Samples')
        self.seed = widgets.IntText(min=1, max=100000000, value=9317915, description='Seed')
        self.btncoradd = widgets.Button(description='+', layout={'width': 'max-content'})
        self.btncorrem = widgets.Button(description='–', layout={'width': 'max-content'})
        self.btnnext = widgets.Button(description='Calculate')
        self.btnlayout = widgets.HBox([widgets.Label('Correlations'),
                                       self.btncoradd, self.btncorrem])
        self.corrs = []

        self.btncoradd.on_click(self.addcor)
        self.btncorrem.on_click(self.remcor)
        super().__init__([self.nsamples, self.seed, self.btnlayout, self.btnnext])

    def addcor(self, change):
        self.corrs.append(Correlation(self.names))
        self.children = [self.nsamples, self.seed, self.btnlayout] + self.corrs + [self.btnnext]

    def remcor(self, change):
        if len(self.corrs) > 0:
            self.corrs.pop(-1)
            self.children = [self.nsamples, self.seed, self.btnlayout] + self.corrs + [self.btnnext]

    def get_corrs(self):
        return [(c.v1.value, c.v2.value, c.val.value) for c in self.corrs]


class OutWidget(widgets.VBox):
    def __init__(self, ucalc):
        self.ucalc = ucalc
        self.view = widgets.Dropdown(value='Summary', options=['Summary', 'Plot', 'Expanded', 'Budget', 'Derivation', 'Convergence'])
        self.view.observe(self.update, names='value')
        plt.ioff()
        self.fig = plt.figure()
        plt.ion()
        self.output = widgets.Output()
        self.update()
        super().__init__([self.view, self.output])

    def update(self, change=None):
        if len(self.ucalc.get_functionnames()) > 1:
            self.view.options = ['Summary', 'Joint Probability', 'Expanded', 'Budget', 'Derivation', 'GUM Validation',
                                 'MC Prob Plot', 'MC Inputs', 'MC Joint Inputs', 'Convergence']
        else:
            self.view.options = ['Summary', 'Expanded', 'Budget', 'Derivation', 'GUM Validation', 
                                 'MC Prob Plot', 'MC Inputs', 'MC Joint Inputs', 'Convergence']

        self.fig.clf()
        with self.output:
            self.output.clear_output()
            if self.view.value == 'Summary':
                display(self.ucalc.out.report(**reportformat))
                self.ucalc.out.plot_pdf(plot=self.fig)
                display(self.fig.canvas)
            elif self.view.value == 'Joint Probability':
                self.ucalc.out.plot_correlation(plot=self.fig, contour=True)
                display(self.fig.canvas)
            elif self.view.value == 'Expanded':
                display(self.ucalc.out.report_expanded(**reportformat))
            elif self.view.value == 'Budget':
                display(self.ucalc.out.report_allinputs(**reportformat))
            elif self.view.value == 'Derivation':
                display(self.ucalc.out.gum.report_derivation(**reportformat))
            elif self.view.value == 'GUM Validation':
                display(self.ucalc.out.report_validity(**reportformat))
            elif self.view.value == 'MC Prob Plot':
                self.ucalc.out.mc.plot_normprob(plot=self.fig)
                display(self.fig.canvas)
            elif self.view.value == 'MC Inputs':
                self.ucalc.out.mc.plot_xhists(plot=self.fig)
                display(self.fig.canvas)
            elif self.view.value == 'MC Joint Inputs':
                self.ucalc.out.mc.plot_xscatter(plot=self.fig)
                display(self.fig.canvas)
            elif self.view.value == 'Convergence':
                self.ucalc.out.mc.plot_converge(plot=self.fig)
                display(self.fig.canvas)

    def message(self, text):
        ''' Show text message (for errors, etc.) '''
        self.output.clear_output()
        with self.output:
            print(text)


class DistributionWidget(widgets.VBox):
    ''' Distribution Edit Widget for Risk Notebooks '''
    def __init__(self):
        self.dist = widgets.Dropdown(options=['normal', 't', 'gamma', 'triangular', 'uniform'], value='normal', description='Distribution')
        self.param1 = widgets.FloatText(value=1, description='Mean')
        self.param2 = widgets.FloatText(value=1, description='Std Dev')
        self.degf = widgets.FloatText(value=10, description='Deg. F.')
        self.degf.layout.visibility = 'hidden'
        self.dist.observe(self.on_distchange)
        super().__init__([self.dist, self.param1, self.param2, self.degf])

    def on_distchange(self, change):
        params = {'normal': ('Mean', 'Std Dev'),
                  't': ('Mean', 'Std Dev'),
                  'gamma': ('a', 'b'),
                  'uniform': ('Mean', 'a'),
                  'triangular': ('Mean', 'a')}.get(self.dist.value)

        self.param1.description = params[0]
        self.param2.description = params[1]
        if self.dist.value == 't':
            self.degf.layout.visibility = 'visible'
        else:
            self.degf.layout.visibility = 'hidden'

    def get_dist(self):
        params = (self.param1.value, self.param2.value)
        if self.dist.value == 'normal':
            dist = stats.norm(*params)
        elif self.dist.value == 't':
            dist = stats.t(df=self.degf.value, loc=params[0], scale=params[1])
        elif self.dist.value == 'gamma':
            dist = stats.gamma(a=params[0], scale=1/params[1])
        elif self.dist.value == 'triangular':
            dist = stats.triang(loc=params[0]-params[1], c=0.5, scale=params[1]*2)
        elif self.dist.value == 'uniform':
            dist = stats.uniform(loc=params[0]-params[1], scale=params[1]*2)
        else:
            raise ValueError
        dist.name = self.dist.value  # For compatibility with suncal.Distribution wrapper
        return dist