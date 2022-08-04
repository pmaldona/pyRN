import { writable } from 'svelte/store';

function createFilePath() {
    const { subscribe, set } = writable("");

    return {
		subscribe,
		open_file: (name) => set(name)
	};
}

export const filename = createFilePath();