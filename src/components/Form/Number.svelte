<script>
    export let bind;
    export let min = undefined;
    export let max = undefined;
    export let step = undefined;
    export let style = "width: 20%;"
</script>

<input type="number" min={min} max={max} step={step} bind:value={bind} style={style}>