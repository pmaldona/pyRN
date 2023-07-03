import { getModal } from '../components/Modal.svelte';

let node_positions = {};

export function openFile(callback) {
    window.electron.open().then(result => {
        let name = result;
        let success = eel.openFile(name)();
        success.then(result => {
            console.log(result);
            if(result == true) {
                callback();
            }
        });
        
    });
}

export function genNetwork() {
    getModal().open();
}

export async function action(random_species, random_vector) {
    console.log(random_species);
    console.log(random_vector);
    let promise = eel.random_network(random_species, random_vector)();
    //await promise.then();
    //let promise_gen = eel.random_network()();
    let randomNetwork = await promise.then(result => {
        return result;
    });
    //getModal().close();
    return true;
}

export function saveNodePositions(nodes) {
    let node_names = Object.keys(nodes);
    node_names.forEach(name => {
        node_positions[name] = [nodes[name].x, nodes[name].y];
    });
}

export function setNodePositions(nodes) {
    let node_names = Object.keys(nodes);
    if(Object.keys(node_positions).length != 0) {
        node_names.forEach(name => {
            if(node_positions[name] == undefined) {
                
            } else {
                nodes[name].x = node_positions[name][0];
                nodes[name].y = node_positions[name][1];
            }
        });
    }
    return nodes;
}