<script>
    export let action;
    export let value;
</script>

<form onsubmit="action(e)">
    <slot></slot>
    <input type="submit" value={value} on:click={action}>
</form>