import tkinter as tk
from tkinter import filedialog,font
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import networkx as nx
from .psec_env import * 
import time 

def get_text_size_in_inches(text_widget):
    # Ensure geometry is updated
    text_widget.update_idletasks()

    # Get widget size in pixels
    width_px = text_widget.winfo_width()
    height_px = text_widget.winfo_height()

    # Get screen DPI (pixels per inch)
    dpi_x = text_widget.winfo_fpixels('1i')
    dpi_y = text_widget.winfo_fpixels('1i')

    # Convert to inches
    width_in = width_px / dpi_x
    height_in = height_px / dpi_y

    return width_in, height_in

def dict_to_networkx(d): 
    G = nx.Graph()

    nodes = list(d.keys()) 
    edges = []
    for k,v in d.items(): 
        for v_ in v: 
            edges.append((k,v_))

    G.add_nodes_from(nodes)
    G.add_edges_from(edges)
    return G 

"""
the main Tkinter application class for xFS user 
interface (xFS UI). 
"""
class PuissecApplication(tk.Frame):

    def __init__(self, master=None):
        tk.Frame.__init__(self, master)

        self.grid()
        self.create_widgets()

        self.se = None 

    def create_widgets(self):
        self.init_primary_window_details() 
        self.set_primary_window_details()

    def init_primary_window_details(self):
        bold_font = font.Font(family="Arial", size=12, weight="bold")
        self.open_button = tk.Button(self, text="OpEn GRaf FiLe", command=self.open_SecEnv)
        self.generate_button = tk.Button(self, text="Generraciones", command=self.generate_SecEnv) 
        self.run_button = tk.Button(self,text="run one", command=self.run_SecEnv_one_time)

    def set_primary_window_details(self): 
        fig = Figure(figsize=(6,6),dpi=100)
        self.canvass = FigureCanvasTkAgg(fig, master=self)
        self.canvass_ax = fig.add_subplot()

        G = nx.Graph()
        nx.draw(G, with_labels=True, font_weight='bold', ax=self.canvass_ax)
        self.canvass.get_tk_widget().grid()
        self.open_button.grid(row=0,column=0) 
        self.generate_button.grid(row=1,column=0)
        self.run_button.grid(row=2,column=0) 

    # TODO: 
    def open_SecEnv(self):
        dirname = filedialog.askdirectory(title="Select a Folder")
        if dirname: 
            raise ValueError("implement this.")
        return 

    def generate_SecEnv(self):
        self.se = default_simple_generate_SecEnv(integer=None,\
            is_naive_hypothesis_type=False,mode_open_info=(1,1))
        self.visualize_SecEnv()
        return 

    def run_SecEnv_one_time(self):
        if type(self.se) == type(None): 
            return 
        
        self.se.run(1.0)
        self.update_SecEnv() 
        return

    def visualize_SecEnv(self): 
        self.canvass_ax.clear() 
        self.G,self.pos = SecEnv.visualize(self.se,ax=self.canvass_ax)  
        self.canvass.draw()  # it needs it in this place
        self.se.preprocess() 

    def update_SecEnv(self):

        g = self.se.sn.d 
        secnodes,cr_loc,ir_loc = self.se.to_vis_data() 
        color_map = default_Puissec_node_color_map(g,secnodes,cr_loc,ir_loc)
        ncs = [cm[1]["color"] for cm in color_map]

        nx.draw_networkx_nodes(self.G, self.pos, node_color = ncs,ax=self.canvass_ax)
        self.canvass.draw() 

def run_puissec_app(): 
    app = PuissecApplication()
    app.master.title('PUISSEC Whendose')  
    app.mainloop()