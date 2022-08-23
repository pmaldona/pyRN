<script>
	import { filename } from './store/FileStore';
	import { network } from './store/NetworkStore';
	import { DataSet } from 'vis-data/peer';
	import Base from './Base.svelte';
  	import Network from './Network.svelte';
	import Hasse from './Hasse.svelte';
	import { hasse } from './store/HasseStore';
    import { synStr } from './store/SynStrStore'; 
	import { protoSyn } from './store/ProtoSynStrStore';
	import Simulation from './Simulation.svelte';

  	const pages = [Base, Network, Hasse, Simulation];

	let page = 0;

	let state = {};

	const unsub_network = network.subscribe();
	const unsub_hasse = hasse.subscribe();
	const unsub_syn_str = synStr.subscribe();
	const unsub_proto_str = protoSyn.subscribe();

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

		return n;
	}

	const generateSynergeticStructure = async () => {
		let promise = eel.gen_synergetic()();
		let s = await promise.then(result => {
			return result;
		});

		return s;
	}

	const generateProtoSynergeticStructure = async () => {
		let promise = eel.gen_protosynergetic()();
		let s = await promise.then(result => {
			return result;
		});

		return s;
	}

	const generateBasicsSpeciesPlot = async () => {
		let promise = eel.plot_basics_sp()();
		let s = await promise.then(result => {
			return result;
		});

		return s;
	}

	const generateBasicsReactionPlot = async () => {
		let promise = eel.plot_basics_r()();
		let s = await promise.then(result => {
			return result;
		});

		return s;
	}

	const generateStoichiometryPlot = async () => {
		let promise = eel.plot_stoichiometry()();
		let s = await promise.then(result => {
			return result;
		});

		return s;
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
	genSynergeticStructure={generateSynergeticStructure}
	genNetwork={generateNetwork}
	genProtoSyn={generateProtoSynergeticStructure}
	setNetwork={setNetwork}
	genBasicsSp={generateBasicsSpeciesPlot}
	genBasicsR={generateBasicsReactionPlot}
	genStoich={generateStoichiometryPlot}
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