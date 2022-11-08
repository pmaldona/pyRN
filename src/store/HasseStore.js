import { writable } from 'svelte/store';
import { generateGraph } from './../actions/drawNetwork';

function createHasse() {
    const { subscribe, set, update } = writable({nodes: [], edges: [], options: {}});
    let hasse = null;

    return {
		subscribe: () => subscribe(h => {
            if(h.nodes) {
                if(h.nodes.length != 0) {
                    hasse = generateGraph('hasse', h);
                }
                
            }
        }),
		set_hasse: (hss) => set(hss),
        get_hasse: () => {
            return hasse;
        },
        update_nodes: (nodes) => update(h => {
            h.nodes = nodes;
        })
	};
}

export const hasse = createHasse();