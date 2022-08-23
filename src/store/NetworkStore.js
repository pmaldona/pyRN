import { writable } from 'svelte/store';
import { generateGraph } from './../actions/drawNetwork';

function createNetwork() {
    const { subscribe, set, update } = writable({nodes: [], edges: [], options: {}, species_count: 0, reaction_count: 0});

    return {
		subscribe: () => subscribe(n => {
            if(n.nodes) {
                if(n.nodes.length != 0) {
                    generateGraph('network', n);
                }
                
            }
            console.log("SUBSCRIBE1")
            return n;
        }),
		set_network: (ntwrk) => set(ntwrk),
        update_nodes: (nodes) => update(n => {
            console.log(nodes);
            console.log(n);
            n.nodes = nodes;
        })
	};
}

export const network = createNetwork();