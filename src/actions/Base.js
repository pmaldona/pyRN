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