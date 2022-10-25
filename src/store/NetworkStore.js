import { writable } from 'svelte/store';
import { generateGraph } from './../actions/drawNetwork';

function createNetwork() {
    const { subscribe, set, update } = writable({nodes: [], edges: [], options: {}, species_count: 0, reaction_count: 0});
    let _net = null;

    return {
		subscribe: () => subscribe(n => {
            if(n.nodes) {
                if(n.nodes.length != 0) {
                    _net = generateGraph('network', n);
                }
                
            }
            console.log("SUBSCRIBE1")
            return n;
        }),
		set_network: (ntwrk) => set(ntwrk),
        get_network: () => {
            return _net;
        },
        update_nodes: (nodes) => update(n => {
            console.log(nodes);
            console.log(n);
            n.nodes = nodes;
        })
	};
}

export const network = createNetwork();