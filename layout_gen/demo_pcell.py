import gdsfactory as gf

@gf.cell
def mzi_with_bend(radius: float=10):
    c = gf.Component()
    #c.add_ref(gf.components.mzi(delta_length=10))
    mzi = c << gf.components.mzi(delta_length=10)
    bend = c << gf.components.bend_circular(radius=radius)
    bend.connect('o1', mzi.ports['o2'])

    c.add_port('o1',port=mzi.ports['o1'])
    c.add_port('o2',port=bend.ports['o2'])
    return c

if __name__ == '__main__':
    c = mzi_with_bend(radius=50)
    c = gf.routing.add_fiber_array(c)
    c.show(show_ports=True)

# c.add_port('o1', port=mzi.ports['o1'])
# c.add_port('o2', port=bend.ports['o2'])
# c.show(show_ports=True)

#gf.read.from_yaml('output.yaml')
