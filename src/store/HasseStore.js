import { writable } from 'svelte/store';
import { generateGraph } from './../actions/drawNetwork';

function createHasse() {
    const { subscribe, set, update } = writable({nodes: [], edges: [], options: {}});

    return {
		subscribe: () => subscribe(h => {
            if(h.nodes) {
                if(h.nodes.length != 0) {
                    generateGraph('hasse', h);
                }
                
            }
        }),
		set_hasse: (hss) => set(hss),
        update_nodes: (nodes) => update(h => {
            h.nodes = nodes;
        })
	};
}

export const hasse = createHasse();