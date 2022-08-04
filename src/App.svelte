<script>
	import { filename } from './store/FileStore';
	import { network } from './store/NetworkStore';
	import { DataSet } from 'vis-data/peer';
	import Base from './Base.svelte';
  	import Network from './Network.svelte';
	import Hasse from './Hasse.svelte';
	import { hasse } from './store/HasseStore';

  	const pages = [Base, Network, Hasse];

	let page = 0;

	let state = {};

	const unsub_network = network.subscribe();
	const unsub_hasse = hasse.subscribe();

	const setName = (name) => {
		filename.open_file(name);
	}

	const openPage = (_p) => {
		page = _p;
		console.log(page);
	}

	const generateHasse = async () => {
		let promise = eel.calculate_orgs()();
		let hasse = await promise.then(result => {
			//console.log(result);
			return result;
		});
		state.hasse = {};
		state.hasse.nodes = new DataSet(hasse.nodes);
		state.hasse.edges = new DataSet(hasse.edges);
		state.hasse.options = new DataSet(hasse.options);
		return hasse;
	}

	const setNetwork = (ntwrk) => {
		console.log(ntwrk)
	}

	const generateNetwork = async () => {
		console.log("genNetwork");
		let promise = eel.gen_network()();
		let n = await promise.then(result => {
			return result;
		});
		//state.network_raw = network;
		// state.network = {};
		// state.network.nodes =new DataSet(network.nodes);
		// state.network.edges = new DataSet(network.edges);
		// state.network.options = new DataSet(network.Options);

		return n;
	}
</script>

<nav>
	<div class="nav-wrapper">
		<ul id="nav-mobile" class="left">
			<li on:click={() => {openPage(0)}}><a>Main</a></li>
			<li on:click={() => {openPage(1)}}><a>Reaction Network</a></li>
			<li on:click={() => {openPage(2)}}><a>Hasse Diagramm</a></li>
			<li on:click={() => {openPage(3)}}><a>Simulation</a></li>
		</ul>
	</div>
</nav>

<svelte:component
  	this={pages[page]}
	{setName}
	genHasse={generateHasse}
	genNetwork={generateNetwork}
	setNetwork={setNetwork}
  	initialValues={state}
/>

<footer>
	<p>{ $filename || "No file opened."}</p>
</footer>

<style>
	footer {
		position: absolute;
		bottom: 0px;
	}
</style>