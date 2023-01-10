import { getModal } from '../components/Modal.svelte';

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

export async function action(with_inflow, random_species, random_reactions, extra, distribution, pr, pp, inflow, outflow, callback) {
    let promise = eel.random_network(with_inflow, random_species, random_reactions, extra, distribution, pr, pp, inflow, outflow)();
    //await promise.then();
    //let promise_gen = eel.random_network()();
    let randomNetwork = await promise.then(result => {
        return result;
    });
    getModal().close();
    if(randomNetwork) {
        callback();
    }
}