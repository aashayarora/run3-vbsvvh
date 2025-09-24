#!/usr/bin/env python
from dataclasses import dataclass, field
import traceback

import ROOT as r
r.EnableImplicitMT()

import matplotlib.pyplot as plt
import mplhep as hep
hep.style.use(hep.style.CMS)

from pathlib import Path

@dataclass
class Hist1D:
    var: str
    xlabel: str
    binning: tuple
    hist_data: list = field(default_factory=list)
    hist_bkg: list = field(default_factory=list)
    hist_sig: list = field(default_factory=list)

@dataclass
class Hist2D:
    xvar: str
    yvar: str
    xlabel: str
    ylabel: str
    xbinning: tuple
    ybinning: tuple
    hist_data: list = field(default_factory=list)
    hist_bkg: list = field(default_factory=list)
    hist_sig: list = field(default_factory=list)

class Plotter:
    def __init__(self, sig=None, bkg=None, data=None, bkg_samples_labels=None, 
                sig_samples_labels=None, cut=None, year=None):
        self.year = year
        self._initialize_dataframes(sig, bkg, data, cut)
        self._setup_sample_labels(bkg_samples_labels, sig_samples_labels)
    
    def _initialize_dataframes(self, sig, bkg, data, cut):
        if data:
            self.df_data = self._create_dataframe(data, cut)
        else:
            self.df_data = None
        
        if sig:
            if isinstance(sig, str):
                self.df_sig = self._create_dataframe(sig, cut)
            elif isinstance(sig, list):
                self.df_sig = [self._create_dataframe(s, cut) for s in sig]
            else:
                self.df_sig = None
        else:
            self.df_sig = None
        
        if bkg:
            self.df_bkg = self._create_dataframe(bkg, cut)
        else:
            self.df_bkg = None
    
    def _create_dataframe(self, source, cut):
        df = r.RDataFrame("Events", source)
        if cut:
            return df.Filter(cut)
        return df
    
    def _setup_sample_labels(self, bkg_samples_labels, sig_samples_labels):
        self.bkg_samples_labels = bkg_samples_labels
        self.sig_samples_labels = sig_samples_labels
        
        if isinstance(self.df_sig, list) and not self.sig_samples_labels:
            raise ValueError("Signal samples labels must be provided for multiple signal samples")
            
        if self.bkg_samples_labels is None and self.df_bkg is not None:
            print("No background labels provided, will use single histogram for background")

    def make_plots(self, hists, density=False, save=True, logy=False, savePath="plots"):
        for histogram in hists:
            self._fill_histogram(histogram)
        
        # Then create all plots
        for hist in hists:
            if isinstance(hist, Hist1D):
                self.plot1D(hist, density=density, logy=logy, save=save)
            elif isinstance(hist, Hist2D):
                self.plot2D(hist, save=save)
    
    def _fill_histogram(self, histogram):
        if isinstance(histogram, Hist1D):
            self._fill_histogram_1d(histogram)
        elif isinstance(histogram, Hist2D):
            self._fill_histogram_2d(histogram)
    
    def _fill_histogram_1d(self, histogram):
        if self.df_data:
            histogram.hist_data = [self.df_data.Histo1D(
                (histogram.var, histogram.var, *histogram.binning), 
                histogram.var, "weight"
            )]
        
        if self.df_sig:
            if isinstance(self.df_sig, list):
                histogram.hist_sig = []
                for sig in self.df_sig:
                    histogram.hist_sig.append(
                        sig.Histo1D((histogram.var, histogram.var, *histogram.binning), 
                                   histogram.var, "weight")
                    )
            else:
                histogram.hist_sig = self.df_sig.Histo1D(
                    (histogram.var, histogram.var, *histogram.binning), 
                    histogram.var, "weight"
                )
        
        if self.df_bkg:
            histogram.hist_bkg = []
            if self.bkg_samples_labels is None:
                histogram.hist_bkg.append(
                    self.df_bkg.Histo1D((histogram.var, histogram.var, *histogram.binning), 
                                       histogram.var, "weight")
                )
            else:
                for sample in self.bkg_samples_labels.keys():
                    histogram.hist_bkg.append(
                        self.df_bkg.Filter(f'type == "{sample}"').Histo1D(
                            (histogram.var, histogram.var, *histogram.binning), 
                            histogram.var, "weight"
                        )
                    )
    
    def _fill_histogram_2d(self, histogram):
        # TODO: Implement 2D histogram filling
        pass

    def plot1D(self, histogram, density=False, save=True, logy=False, savePath="plots"):
        try:
            hist_ratio = self._create_ratio_histogram(histogram)
            fig, ax_main, ax_ratio = self._setup_figure_axes(hist_ratio is not None)
            self._plot_background_histograms(histogram, ax_main, density)
            self._plot_data_histograms(histogram, ax_main, density)
            self._plot_signal_histograms(histogram, ax_main, density)

            if logy:
                ax_main.set_yscale("log")
                ax_main.set_ylim(0.1, None)
            
            if hist_ratio:
                hep.histplot(hist_ratio, color="black", ax=ax_ratio, histtype="errorbar", density=density)
            
            year = self.year if self.year else "Run3"
            hep.cms.label("Preliminary", data=True, year=year, ax=ax_main)
            
            self._configure_axes(ax_main, ax_ratio, histogram, hist_ratio is not None)
            if save:
                self._save_plot(fig, histogram.var, savePath)
                
            plt.close(fig)
            
        except Exception as e:
            print(f"Error plotting {histogram.var}: {e}")
            traceback.print_exc()

    def _create_ratio_histogram(self, histogram):
        if (self.df_data and self.df_bkg and 
            histogram.hist_data and histogram.hist_bkg):
            hist_ratio = histogram.hist_data[0].GetValue().Clone()
            hist_bkg_total = histogram.hist_bkg[0].GetValue().Clone()
            for hist in histogram.hist_bkg[1:]:
                hist_bkg_total.Add(hist.GetValue())
            hist_ratio.Divide(hist_bkg_total)
            return hist_ratio
        return None

    def _setup_figure_axes(self, has_ratio):
        if has_ratio:
            fig, ax = plt.subplots(2, 1, gridspec_kw={"height_ratios": (4, 1)})
            return fig, ax[0], ax[1]
        else:
            fig, ax = plt.subplots(1, 1)
            return fig, ax, None

    def _plot_signal_histograms(self, histogram, ax, density):
        if not histogram.hist_sig:
            return
            
        if isinstance(histogram.hist_sig, list):
            hist_values = [h.GetValue() for h in histogram.hist_sig]
            for hist in hist_values:
                hist.Scale(1000)  # Scale by 1000 for visibility
                
            for i, hist in enumerate(hist_values):
                hep.histplot(
                    hist, ax=ax, histtype="step", 
                    label=f"Signal {self.sig_samples_labels[i]} x 1000", 
                    linewidth=3, yerr=False, density=density
                )
        else:
            # Single signal case
            hist_value = histogram.hist_sig.GetValue()
            hist_value.Scale(1000)
            hep.histplot(
                hist_value, ax=ax, histtype="step", 
                label=f"Signal x 1000", linewidth=3, yerr=False, 
                density=density, color="red"
            )

    def _plot_background_histograms(self, histogram, ax, density):
        if not histogram.hist_bkg:
            return
            
        if self.bkg_samples_labels is None:
            # Single background case
            hep.histplot(
                [h.GetValue() for h in histogram.hist_bkg], 
                ax=ax, histtype="fill", label="Background", 
                density=density
            )
        else:
            # Multiple background samples case
            hep.histplot(
                [h.GetValue() for h in histogram.hist_bkg], 
                ax=ax, histtype="fill", stack=True, 
                label=list(self.bkg_samples_labels.values()), 
                density=density
            )

    def _plot_data_histograms(self, histogram, ax, density):
        if histogram.hist_data:
            hep.histplot(
                [h.GetValue() for h in histogram.hist_data], 
                label="Data", ax=ax, histtype="errorbar", 
                color="black", density=density
            )

    def _configure_axes(self, ax_main, ax_ratio, histogram, has_ratio):
        ax_main.legend()
        ax_main.set_ylabel("Events")
        
        if has_ratio:
            ax_main.set_xlabel("")
            ax_ratio.set_xlabel(histogram.xlabel)
            ax_ratio.set_ylabel("Data / MC")
            ax_ratio.set_ylim(0.8, 1.2)
            ax_ratio.axhline(1, color="black", linestyle="--")
        else:
            ax_main.set_xlabel(histogram.xlabel)

    def _save_plot(self, fig, plot_name, savePath):
        Path(savePath).mkdir(parents=True, exist_ok=True)
        fig.tight_layout()
        plot_path = f"{savePath}/{plot_name}.png"
        plt.savefig(plot_path)
        print(f"Saved plot to {plot_path}")

    def plot2D(self, histogram, save=True, savePath="plots"):
        # TODO: Implement the 2D plotting functionality
        raise NotImplementedError("2D plotting not yet implemented")