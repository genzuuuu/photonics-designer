# Layout Automated Generator
This module is aimed to achieve the function of layout automated generator in the whole workflow.

The first necessary accessibility is the utilization of [Klayout](https://www.klayout.de/) which is an integrated circuit layout editor in photonics design domain.
*demo_pcell.py* is the file which converts the programming language to the layout .gds file via [gdsfactory package](https://gdsfactory.github.io/gdsfactory/). This demo file gives the example of defining the canvas, adding different components and adding ports.

*20240207 llama_index_KG.ipynb* and *20240214 layout_gen with llama_vector_store* gives the procedure to creating specific-domain knowledge graph by importing related information. Based on this graph, some tests are implenmented as well to generate the code following the function in gdsfactory.

The combination of layout and simulation requires the file conversion between Klayout and [Lumerical INTERCONNECT](https://www.ansys.com/products/optics/interconnect). The script command *exportnetlist(filename)* can export the netlist file in spice format. *20240220 spi_to_yaml* is possible to convert the file type to .yaml format. Since the lanuage rule is different, this conversion results in the .yaml file is unreadable under gdsfactory package supporting which requires further debugging.
